package jq255s

import (
	"math/bits"
	"github.com/doubleodd/go-jq255/internal/scalar"
)

// This file defines types and functions for jq255s scalars, i.e.
// integers modulo r = 2^254 + 56904135270672826811114353017034461895
//
// Unless explicitly documented, all functions here are constant-time.

// Jq255sScalar is the type for an integer modulo the prime order of the
// jq255s group. Default value is zero.
type Jq255sScalar [4]uint64

// Decode a scalar from exactly 32 bytes. Returned value is:
//   1   scalar properly decoded, value is not zero
//   0   scalar properly decoded, value is zero
//  -1   source bytes were not a valid scalar encoding
// If the decoding fails, then the scalar value is forced to zero.
func (s *Jq255sScalar) Decode(src []byte) int {
	return scalar.Decode((*[4]uint64)(s), src, &jq255sOrder)
}

// Decode a scalar from some bytes. All provided bytes are read and
// interpreted as an integer in unsigned little endian convention, which
// is reduced modulo the curve subgroup order. This process cannot fail.
func (s *Jq255sScalar) DecodeReduce(src []byte) {
	scalar.DecodeReduce((*[4]uint64)(s), src, jq255sModrReduce384Partial)
}

// Encode a scalar into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (s *Jq255sScalar) Encode(dst []byte) []byte {
	return scalar.Encode(dst, (*[4]uint64)(s), jq255sModrReduce256)
}

// Encode a scalar into exactly 32 bytes; the newly allocated slice
// is returned.
func (s *Jq255sScalar) Bytes() [32]byte {
	return scalar.ToBytes((*[4]uint64)(s), jq255sModrReduce256)
}

// Compare a scalar with zero. Returned value is 1 if the scalar is zero,
// 0 otherwise.
func (s *Jq255sScalar) IsZero() int {
	var t [4]uint64
	jq255sModrReduce256(&t, (*[4]uint64)(s))
	z := t[0] | t[1] | t[2] | t[3]
	return int(1 - ((z | -z) >> 63))
}

// Compare two scalars together. Returned value is 1 if the scalars are
// equal to each other, 0 otherwise.
func (s *Jq255sScalar) Equal(a *Jq255sScalar) int {
	var t Jq255sScalar
	t.Sub(s, a)
	return t.IsZero()
}

// Scalar addition: s is set to a + b (mod r).
// A pointer to s is returned.
func (s *Jq255sScalar) Add(a, b *Jq255sScalar) *Jq255sScalar {
	scalar.Add((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), jq255sModrReduce256Partial)
	return s
}

// Scalar subtraction: s is set to a - b (mod r).
// A pointer to s is returned.
func (s *Jq255sScalar) Sub(a, b *Jq255sScalar) *Jq255sScalar {
	scalar.Sub((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), jq255sModrReduce256Partial, &jq255sOrder)
	return s
}

// Scalar negation: s is set to -a (mod r).
// A pointer to s is returned.
func (s *Jq255sScalar) Neg(a *Jq255sScalar) *Jq255sScalar {
	var t = [4]uint64 { 0, 0, 0, 0 }
	scalar.Sub((*[4]uint64)(s), &t, (*[4]uint64)(a), jq255sModrReduce256Partial, &jq255sOrder)
	return s
}

// Scalar multiplication: s is set to a*b (mod r).
// A pointer to s is returned.
func (s *Jq255sScalar) Mul(a, b *Jq255sScalar) *Jq255sScalar {
	scalar.Mul((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), jq255sModrReduce384Partial)
	return s
}

// Group order is r = 2^254 + r0, with:
//    r0 = 56904135270672826811114353017034461895
const jq255s_r0_lo uint64 = 0xDCF2AC65396152C7
const jq255s_r0_hi uint64 = 0x2ACF567A912B7F03
const jq255s_r_top uint64 = 0x4000000000000000

const jq255s_r0x4_lo uint64 = 0x73CAB194E5854B1C
const jq255s_r0x4_hi uint64 = 0xAB3D59EA44ADFC0F

var jq255sOrder = [4]uint64 { jq255s_r0_lo, jq255s_r0_hi, 0, jq255s_r_top }

var jq255s_r2 = [8]uint64 {
	0xA31F34E2739216B1, 0x86A297C9835B5211,
	0x95DCE66BF04303AD, 0x8728B04D2F0F9E3C,
	0xEE7956329CB0A963, 0x1567AB3D4895BF81,
	0x0000000000000000, 0x1000000000000000,
}

// Given input 'a' (up to 2^286-1), perform a partial reduction modulo r;
// output (into 'd') fits on 255 bits and is lower than 2*r. The
// high bits of 'a' are provided as extra parameter ah.
func jq255sModrReduce256PartialWithExtra(d, a *[4]uint64, ah uint64) {
	// Truncate to 254 bits and get extra bits into ah.
	ah = (ah << 2) | (a[3] >> 62)

	// Compute ah*r0 into u0:u1:u2.
	u1, u0 := bits.Mul64(ah, jq255s_r0_lo)
	u2, lo := bits.Mul64(ah, jq255s_r0_hi)
	var cc uint64
	u1, cc = bits.Add64(u1, lo, 0)
	u2 += cc

	// 2^254 = -r0 mod r
	d[0], cc = bits.Sub64(a[0], u0, 0)
	d[1], cc = bits.Sub64(a[1], u1, cc)
	d[2], cc = bits.Sub64(a[2], u2, cc)
	d[3], cc = bits.Sub64(a[3] & 0x3FFFFFFFFFFFFFFF, 0, cc)

	// If we got a borrow, then we must add back r. Since ah*r0 < 2^192,
	// the result will be nonnegative, but less than r.
	m := -cc
	d[0], cc = bits.Add64(d[0], m & jq255sOrder[0], 0)
	d[1], cc = bits.Add64(d[1], m & jq255sOrder[1], cc)
	d[2], cc = bits.Add64(d[2], m & jq255sOrder[2], cc)
	d[3] = d[3] + (m & jq255sOrder[3]) + cc
}

// Partial reduction ensures that the output is fits on 255 bits and is
// less than 2*r.
func jq255sModrReduce256Partial(d, a *[4]uint64) {
	jq255sModrReduce256PartialWithExtra(d, a, 0)
}

// Given a partially reduced input 'a' (less than 2*r), finish reduction
// (conditional subtraction of r).
func jq255sModrReduce256Finish(d, a *[4]uint64) {
	// Try to subtract r.
	var t [4]uint64
	var cc uint64
	t[0], cc = bits.Sub64(a[0], jq255s_r0_lo, 0)
	t[1], cc = bits.Sub64(a[1], jq255s_r0_hi, cc)
	t[2], cc = bits.Sub64(a[2], 0, cc)
	t[3], cc = bits.Sub64(a[3], jq255s_r_top, cc)

	// If the result is nonnegative, then keep it; otherwise, use the
	// original value.
	m := -cc
	for i := 0; i < 4; i ++ {
		d[i] = t[i] ^ (m & (a[i] ^ t[i]))
	}
}

// Perform full reduction of a scalar.
func jq255sModrReduce256(d, a *[4]uint64) {
	jq255sModrReduce256Partial(d, a)
	jq255sModrReduce256Finish(d, d)
}

// Given a 384-bit input 'a', perform a partial reduction modulo r;
// output fits on 255 bits and is less than 2*r.
func jq255sModrReduce384Partial(d *[4]uint64, a *[6]uint64) {
	// Multiply the high third (a4:a5) by 4*r0 into tw.
	var t1, t2 [2]uint64
	var tw [4]uint64
	t1[0] = jq255s_r0x4_lo
	t1[1] = jq255s_r0x4_hi
	t2[0] = a[4]
	t2[1] = a[5]
	scalar.Mul128x128(&tw, &t1, &t2)

	// Subtract 4*r0*ah from the low part of 'a', then
	// add back 4*r. Since 4*r0 =~ 2^127.42, the result may be
	// slightly above 2^257, but will fit on 258 bits.
	var cc uint64
	tw[0], cc = bits.Sub64(a[0], tw[0], 0)
	tw[1], cc = bits.Sub64(a[1], tw[1], cc)
	tw[2], cc = bits.Sub64(a[2], tw[2], cc)
	tw[3], cc = bits.Sub64(a[3], tw[3], cc)
	tw4 := -cc
	tw[0], cc = bits.Add64(tw[0], jq255s_r0x4_lo, 0)
	tw[1], cc = bits.Add64(tw[1], jq255s_r0x4_hi, cc)
	tw[2], cc = bits.Add64(tw[2], 0, cc)
	tw[3], cc = bits.Add64(tw[3], 0, cc)
	tw4 += (cc + 1)

	// Perform partial reduction.
	jq255sModrReduce256PartialWithExtra(d, &tw, tw4)
}

// Recode a scalar with 5-bit Booth encoding. Output is a sequence of
// small integers in the -15..+16 range such that:
//   a = \sum_{i=0}^{51} d[i]*2^(5*i)
// Top digit d[51] is nonnegative. If the input value is less than 2^255,
// then the top digit can only be 0 or 1.
// Each output digit is encoded in a byte as sign+mantissa: the low 5 bits
// of the byte are the absolute value of the digit (in the 0..16 range),
// and the high bit of the byte is set to 1 for a negative digit, 0 otherwise.
// When the digit is 0, the function may encode it as -0 (0x80) or +0 (0x00)
// (the top digit d[51] cannot be -0, only +0).
func (a *Jq255sScalar) recode5(d *[52]byte) {
	scalar.Recode5(d, (*[4]uint64)(a))
}
