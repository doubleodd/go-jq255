package jq255e

import (
	"bytes"
	"errors"
	"io"
	"encoding/binary"
	"crypto"
	cryptorand "crypto/rand"
	"golang.org/x/crypto/blake2s"
)

// This file implements high-level operations over jq255e:
//
//   - Key pair generation
//   - Key exchange (ECDH)
//   - Signature generation and verification
//   - Hash-to-curve

// A private key structure contains a private key, i.e. a non-zero
// scalar for jq255e. For efficiency reasons, it internally caches a
// copy of the public key as well.
type Jq255ePrivateKey struct {
	d Jq255eScalar
	pub Jq255ePoint
	epub [32]byte
}

// A public key structure contains a non-neutral group element.
type Jq255ePublicKey struct {
	pub Jq255ePoint
	epub [32]byte
}

// Test whether a public key is equal to another.
func (pk Jq255ePublicKey) Equal(other crypto.PublicKey) bool {
	pk2, ok := other.(Jq255ePublicKey)
	if !ok {
		return false
	}
	var t byte = 0
	for i := 0; i < 32; i ++ {
		t |= pk.epub[i] ^ pk2.epub[i]
	}
	return t == 0
}

// Decode a private key from bytes. This function expects exactly
// 32 bytes. If the provided slice does not have length exactly 32,
// or if what it contains is not the canonical encoding of a valid
// non-zero scalar for jq255e, then this function returns nil and an
// error.
func Jq255eDecodePrivateKey(src []byte) (*Jq255ePrivateKey, error) {
	if len(src) != 32 {
		return nil, errors.New("Invalid private key")
	}
	sk := new(Jq255ePrivateKey)
	if sk.d.Decode(src) != 1 {
		return nil, errors.New("Invalid private key")
	}
	sk.pub.MulGen(&sk.d)
	sk.pub.Encode(sk.epub[:0])
	return sk, nil
}

// Encode a private key into bytes. The private key (exactly 32 bytes)
// is appended to the provided slice. If 'dst' has enough capacity, then
// it is returned; otherwise, a new slice is allocated, and receives
// the concatenation of the current contents of 'dst' and the encoded
// private key.
func (sk *Jq255ePrivateKey) Encode(dst []byte) []byte {
	return sk.d.Encode(dst)
}

// Get the public key corresponding to a given private key.
func (sk *Jq255ePrivateKey) Public() *Jq255ePublicKey {
	pk := new(Jq255ePublicKey)
	pk.pub.Set(&sk.pub)
	copy(pk.epub[:], sk.epub[:])
	return pk
}

// Decode a public key from bytes. This function expects exactly
// 32 bytes. If the provided slice does not have length exactly 32,
// or if what it contains is not the canonical encoding of a valid
// non-neutral jq255e element, then this function returns nil and an
// error.
func Jq255eDecodePublicKey(src []byte) (*Jq255ePublicKey, error) {
	if len(src) != 32 {
		return nil, errors.New("Invalid public key")
	}
	pk := new(Jq255ePublicKey)
	if pk.pub.Decode(src) != 1 {
		return nil, errors.New("Invalid public key")
	}
	copy(pk.epub[:], src)
	return pk, nil
}

// Encode a public key into bytes. The public key (exactly 32 bytes)
// is appended to the provided slice. If 'dst' has enough capacity, then
// it is returned; otherwise, a new slice is allocated, and receives
// the concatenation of the current contents of 'dst' and the encoded
// public key.
func (pk *Jq255ePublicKey) Encode(dst []byte) []byte {
	n := len(dst)
	n2 := n + 32
	var dst2 []byte
	if cap(dst) >= n2 {
		dst2 = dst[:n2]
	} else {
		dst2 = make([]byte, n2)
		copy(dst2, dst)
	}
	copy(dst2[n:], pk.epub[:])
	return dst2
}

// Key pair generation with jq255e: from a random source 'rand', a
// private key (a scalar value) and the corresponding public key (group
// element) are generated. The random source MUST be cryptographically
// secure. If 'rand' is nil, then crypto/rand.Reader is used (this is
// the recommended way).
func Jq255eGenerateKeyPair(rand io.Reader) (*Jq255ePrivateKey, error) {
	// We obtain 32 bytes from the random source and decode them
	// as a scalar with modular reduction. The group order is very
	// close to 2^254 (difference is less than 2^127) so the bias
	// is negligible.
	if rand == nil {
		rand = cryptorand.Reader
	}
	var bb [32]byte
	if _, err := io.ReadFull(rand, bb[:]); err != nil {
		return nil, err
	}
	sk := new(Jq255ePrivateKey)
	sk.d.DecodeReduce(bb[:])

	// In the utterly improbable case that we got a zero here, we
	// replace it with 1.
	var bb2 [32]byte
	bb2[0] = byte(sk.d.IsZero())
	var d2 Jq255eScalar
	d2.Decode(bb2[:])
	sk.d.Add(&sk.d, &d2)

	// Compute public key.
	sk.pub.MulGen(&sk.d)
	sk.pub.Encode(sk.epub[:0])

	return sk, nil
}

// Key exchange with jq255e: given our private key, and the public key
// from the peer, a 32-byte shared secret is produced. The peer's public
// key is provided encoded; it should have length exactly 32 bytes. If
// the provided sequence of bytes has not length exactly 32 bytes, or if
// it is not otherwise a valid jq255e point encoding, then the key
// exchange fails. On failure, a byte sequence of the requested length
// is still generated; that byte sequence is not predictable by
// outsiders, and cannot be distinguished from the output of a
// successful ECDH exchange by outsiders. This is meant to support rare
// protocols where exchanged keys are not public, and the exchange
// should not have any validation semantics. The 'ok' returned value has
// value 1 on success, 0 on error (an 'int' is used to promote
// constant-time processing).
func Jq255eKeyExchange(sk *Jq255ePrivateKey, peer_pk []byte) (secret [32]byte, ok int) {
	// Decode peer point; remember if it worked.
	var P Jq255ePoint
	var venc int
	var eppk [32]byte
	if len(peer_pk) == 32 {
		venc = P.Decode(peer_pk)
		copy(eppk[:], peer_pk)
	} else {
		P.Neutral()
		venc = -1
	}
	// Input point is valid if and only if it decoded properly AND it
	// was not the neutral (i.e. venc == 1, not 0 or -1).
	ok = int(-uint32(venc) >> 31)

	// ECDH
	var P2 Jq255ePoint
	P2.Mul(&P, &sk.d)

	// Generate secret as BLAKE2s-256 over the concatenation of:
	//  - our public key and the peer public key, both encoded;
	//    first one is the numerically lowest, when values are
	//    interpreted as integers in unsigned big-endian convention
	//    (if the peer public key does not have length 32 bytes
	//    exactly, then our public key goes first)
	//  - a byte of value 0x53 (success) or 0x46 (failure)
	//  - the encoded shared point (on success) or our private key
	//    (on failure)
	sh, _ := blake2s.New256(nil)

	if len(peer_pk) == 32 {
		// Get the two public keys in numerical order. Our public
		// key is in sk.epub; the peer public key is in eppk.
		var cc uint = 0
		for i := 31; i >= 0; i -- {
			x := uint(sk.epub[i]) - uint(eppk[i]) - cc
			cc = (x >> 8) & 1
		}
		// If cc == 1, then sk.epub < eppk; otherwise, cc == 0 and
		// sk.epub >= eppk.
		var bb2 [64]byte
		m := byte(-cc)
		for i := 0; i < 32; i ++ {
			bb2[i] = (sk.epub[i] & m) | (eppk[i] & ^m)
			bb2[32 + i] = (sk.epub[i] & ^m) | (eppk[i] & m)
		}
		sh.Write(bb2[:])
	} else {
		// If the peer public key has not length exactly 32 bytes
		// then our public key is hashed first.
		sh.Write(sk.epub[:])
		sh.Write(peer_pk)
	}

	// Get bb[0] == 0x53 on success, 0x46 otherwise.
	var bb [33]byte
	var esk [32]byte
	okm := -byte(ok)
	bb[0] = 0x46 + (okm & 0x0D)
	// bb[1..32] is either the shared point, or our own private key.
	P2.Encode(bb[1:1])
	sk.Encode(esk[:0])
	for i := 0; i < 32; i ++ {
		bb[1 + i] ^= ^okm & (bb[1 + i] ^ esk[i])
	}
	sh.Write(bb[:])

	// Generate and return secret.
	sh.Sum(secret[:0])
	return
}

// Get the label string for the provided hash identifier.
func getHashLabel(opts crypto.SignerOpts) (label []byte, err error) {
	switch opts.HashFunc() {
	case crypto.Hash(0):
		label = []byte("")
	case crypto.SHA224:
		label = []byte("sha224")
	case crypto.SHA256:
		label = []byte("sha256")
	case crypto.SHA384:
		label = []byte("sha384")
	case crypto.SHA512:
		label = []byte("sha512")
	case crypto.SHA512_224:
		label = []byte("sha512224")
	case crypto.SHA512_256:
		label = []byte("sha512256")
	case crypto.SHA3_224:
		label = []byte("sha3224")
	case crypto.SHA3_256:
		label = []byte("sha3256")
	case crypto.SHA3_384:
		label = []byte("sha3384")
	case crypto.SHA3_512:
		label = []byte("sha3512")
	case crypto.BLAKE2s_256:
		label = []byte("blake2s")
	case crypto.BLAKE2b_512:
		label = []byte("blake2b")
	// TODO: add label for BLAKE3, then the constant is defined in
	// Go's standard library.
	default:
		return nil, errors.New("Unknown hash identifier")
	}
	return label, nil
}

// Compute the "challenge" in a signature process.
func mkChallenge(Renc []byte, epub []byte, hashLabel []byte, data []byte) [16]byte {
	sh, _ := blake2s.New256(nil)
	sh.Write(Renc)
	sh.Write(epub)
	if len(hashLabel) == 0 {
		sh.Write([]byte{0x52})
	} else {
		sh.Write([]byte{0x48})
		sh.Write(hashLabel)
		sh.Write([]byte{0x00})
	}
	sh.Write(data)
	var tmp [32]byte
	var challenge [16]byte
	sh.Sum(tmp[:0])
	copy(challenge[:], tmp[:16])
	return challenge
}

// Schnorr signature with jq255e. The data to sign ('data') may be either
// raw data, or a hash value. The 'opts' parameter specifies the hash
// function that was used to pre-hash the data (use crypto.Hash(0) for
// raw data).
//
// The signature process is deterministic, for a given 'seed' value. The
// seed may be nil, which is equivalent to a seed of length 0. How the
// seed contents were chosen does not impact the algorithmic security of
// the signature (using a non-repeating seed may increase implementation
// robustness against some classes of physical attacks).
//
// The signature is returned as a newly allocated slice. Its length is
// exactly 48 bytes. An error is reported if the hash function
// identified by 'opts' is not known.
func (sk *Jq255ePrivateKey) signWithSeed(seed []byte, data []byte, opts crypto.SignerOpts) (signature []byte, err error) {
	// Get hash label.
	var hashLabel []byte
	hashLabel, err = getHashLabel(opts)
	if err != nil {
		return nil, err
	}

	// Generate the per-signature k value.
	sh, _ := blake2s.New256(nil)
	var bb [32]byte
	sk.d.Encode(bb[:0])
	sh.Write(bb[:32])
	sh.Write(sk.epub[:])
	var tt [8]byte
	if seed == nil {
		sh.Write(tt[:8])
	} else {
		binary.LittleEndian.PutUint64(tt[:], uint64(len(seed)))
		sh.Write(tt[:8])
		sh.Write(seed)
	}
	if len(hashLabel) == 0 {
		sh.Write([]byte{0x52})
	} else {
		sh.Write([]byte{0x48})
		sh.Write(hashLabel)
		sh.Write([]byte{0x00})
	}
	sh.Write(data)
	sh.Sum(bb[:0])
	var k Jq255eScalar
	k.DecodeReduce(bb[:])

	// R = k*G
	var R Jq255ePoint
	R.MulGen(&k)
	var R_enc [32]byte
	R.Encode(R_enc[:0])

	// Compute challenge c.
	c := mkChallenge(R_enc[:], sk.epub[:], hashLabel, data)

	// s = k + c*d
	var s Jq255eScalar
	s.DecodeReduce(c[:])
	s.Mul(&s, &sk.d).Add(&s, &k)

	// Signature is (c,s)
	signature = make([]byte, 48)
	copy(signature[:16], c[:])
	s.Encode(signature[:16])
	return
}

// Schnorr signature with jq255e. The data to sign ('data') may be either
// raw data, or a hash value. The 'opts' parameter specifies the hash
// function that was used to pre-hash the data (use crypto.Hash(0) for
// raw data).
//
// If 'rand' is nil, then the signature is deterministic (this is safe).
// If 'rand' is not nil, then 32 bytes are read from it and used to
// complement the internal per-signature nonce generation process,
// making the signature non-deterministic, in case a specific protocol
// requires this property. Non-deterministic signatures might also
// improve implementation robustness against some kinds of physical
// attacks (in particular fault attacks). It is not necessary that the
// extra randomness returned by 'rand' has high quality; the security of
// the signature will be maintained in all case, even if that data is
// fully predictable.
//
// The signature is returned as a newly allocated slice. Its length is
// exactly 48 bytes. An error is reported if 'rand' is not nil but a
// read attempt returns an error. An error is also reported if the hash
// function identified by 'opts' is not known.
func (sk *Jq255ePrivateKey) Sign(rand io.Reader, data []byte, opts crypto.SignerOpts) (signature []byte, err error) {
	if rand == nil {
		return sk.signWithSeed(nil, data, opts)
	} else {
		var seed [32]byte
		_, err = io.ReadFull(rand, seed[:])
		if err != nil {
			return nil, err
		}
		return sk.signWithSeed(seed[:], data, opts)
	}
}

// Verify a signature on a message, relatively to a public key.
//
// The message data is provided in 'data'. This is interpreted as raw data
// if opts is crypto.Hash(0)); otherwise, it will be considered to be
// pre-hashed data, processed with the hash function identified by opts.
//
// Returned value is true if the hash function is recognized, and the
// signature is valid relatively to the provided public key. In all other
// cases, false is returned.
//
// This function is not constant-time, under the assumption that public
// keys and signatures are public data.
func (pk *Jq255ePublicKey) VerifyVartime(data []byte, opts crypto.SignerOpts, sig []byte) bool {
	// Valid signatures have length 48 bytes exactly.
	if len(sig) != 48 {
		return false
	}

	// Sanity check: public keys cannot be the neutral point.
	if pk.pub.IsNeutral() != 0 {
		return false
	}

	// Get the hash function label string. If unrecognized, then we
	// have to reject the signature.
	hashLabel, _ := getHashLabel(opts)
	if hashLabel == nil {
		return false
	}

	// The signature splits into the challenge and a scalar. Decode
	// the scalar; reject the signature if that value is not canonical.
	var c [2]uint64
	c[0] = binary.LittleEndian.Uint64(sig[:8])
	c[1] = binary.LittleEndian.Uint64(sig[8:])
	var s Jq255eScalar
	if s.Decode(sig[16:]) < 0 {
		return false
	}

	// Compute R = s*G - c*pk.
	var Q Jq255ePoint
	Q.Neg(&pk.pub)
	var R Jq255ePoint
	R.Mul128AddMulGenVartime(&Q, &c, &s)

	// Recompute the challenge and compare it with the received value.
	var R_enc [32]byte
	R.Encode(R_enc[:0])
	c2 := mkChallenge(R_enc[:], pk.epub[:], hashLabel, data)
	return bytes.Equal(sig[:16], c2[:])
}

// Hash some input data into a curve point. The input data ('data') is
// either raw or pre-hashed, as identified by the opts parameter. If
// the hash function identifier is not recognized, then an error is
// returned. Otherwise, the output point is returned.
func Jq255eHashToCurve(data []byte, opts crypto.SignerOpts) (*Jq255ePoint, error) {
	// Get internal identifier for the hash function.
	label, err := getHashLabel(opts)
	if err != nil {
		return nil, err
	}

	// Hash the input into two 32-byte blobs with BLAKE2s.
	sh, _ := blake2s.New256(nil)
	var bb1 [32]byte
	var bb2 [32]byte
	if len(label) == 0 {
		sh.Write([]byte{0x01, 0x52})
		sh.Write(data)
		sh.Sum(bb1[:0])
		sh.Reset()
		sh.Write([]byte{0x02, 0x52})
		sh.Write(data)
		sh.Sum(bb2[:0])
	} else {
		sh.Write([]byte{0x01, 0x48})
		sh.Write(label)
		sh.Write([]byte{0x00})
		sh.Write(data)
		sh.Sum(bb1[:0])
		sh.Reset()
		sh.Write([]byte{0x02, 0x48})
		sh.Write(label)
		sh.Write([]byte{0x00})
		sh.Write(data)
		sh.Sum(bb2[:0])
	}

	// Map each blob into a curve point, and add them together.
	var P1, P2 Jq255ePoint
	P1.MapBytes(bb1[:])
	P2.MapBytes(bb2[:])
	return NewJq255ePoint().Add(&P1, &P2), nil
}
