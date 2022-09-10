package jq255e

import (
	gf "github.com/doubleodd/go-jq255/internal/field"
	"github.com/doubleodd/go-jq255/internal/scalar"
)

// This file implements operations on curve points for jq255e
// (specifically, on elements of the prime order group defined over
// jq255e).
//
// API: a point is represented in memory by a Jq255ePoint structure.
// Contents of such a structure are opaque. These structures are
// mutable; the various functions such as setAdd() modify the point
// on which they are called. It is always acceptable to also use the
// destination structure as one of the operands. All such functions
// return a pointer to the structure on which they were called, so
// that calls may be syntactically chained.
//
// A point can be encoded to, and decoded from, a sequence of 32 bytes.
// Encoding is unique and verified. Decoding an invalid sequence of
// bytes yields an error. The group neutral element can be encoded (as
// a sequence of 32 bytes of value 0x00). Whether the neutral element
// is acceptable or not when decoding depends on usage context; the
// decoding functions returns a specific flag value when the point is
// the neutral, that the caller may test explicitly.
//
// Unless explicitly documented, all functions here are constant-time.

// Jq255ePoint is the type for a jq255e point.
//
// Default value for a point structure is not valid. The NewJq255ePoint()
// function makes sure to return only initialized structures. If allocating
// a point structure manually, make sure to properly set it to a valid point
// before using it as source.
type Jq255ePoint struct {
	// Internally, we use extended (e,u) coordinates, which have
	// complete and efficient formulas.
	e, z, u, t gf.GF255e
}

// Preallocated neutral point. Do not modify.
var jq255eNeutral = Jq255ePoint {
	e: gf.GF255e { 1, 0, 0, 0 },
	z: gf.GF255e { 1, 0, 0, 0 },
	u: gf.GF255e { 0, 0, 0, 0 },
	t: gf.GF255e { 0, 0, 0, 0 },
}

// Preallocated conventional generator point. Do not modify.
var jq255eGenerator = Jq255ePoint {
	e: gf.GF255e { 3, 0, 0, 0 },
	z: gf.GF255e { 1, 0, 0, 0 },
	u: gf.GF255e { 1, 0, 0, 0 },
	t: gf.GF255e { 1, 0, 0, 0 },
}

// Create a new point. The point is set to the group neutral element (N).
func NewJq255ePoint() *Jq255ePoint {
	P := new(Jq255ePoint)
	*P = jq255eNeutral
	return P
}

// Set the point P to the neutral element (N).
// A pointer to this structure is returned.
func (P *Jq255ePoint) Neutral() *Jq255ePoint {
	*P = jq255eNeutral
	return P
}

// Set the point P to the conventional generator (G).
// A pointer to this structure is returned.
func (P *Jq255ePoint) Generator() *Jq255ePoint {
	*P = jq255eGenerator
	return P
}

// Encode a point into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (P *Jq255ePoint) Encode(dst []byte) []byte {
	// Encoded value is u or -u: if e is "negative" (least significant
	// bit of the encoding of e is 1), then we use -u, otherwise
	// we use u.
	var iz gf.GF255e
	iz.Inv(&P.z)
	var r gf.GF255e
	r.Mul(&iz, &P.e)
	nn := r.IsNegative()
	r.Mul(&iz, &P.u)
	r.CondNeg(&r, nn)
	return r.Encode(dst)
}

// Encode a point into exactly 32 bytes.
func (P *Jq255ePoint) Bytes() [32]byte {
	var d [32]byte
	P.Encode(d[:0])
	return d
}

// Test whether a given chunk of 32 bytes is a valid representation of
// a jq255e group element. This is faster than actually decoding it.
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral group element
//   -1   invalid encoding
func Jq255eCheckPoint(src []byte) int {
	// The bytes can be decoded as a valid point if and only if all of
	// the following hold:
	//  - The bytes can be decoded as a field element u.
	//  - The value (a^2-4*b)*u^4 - 2*a*u^2 + 1 is a square.
	// On jq255e, a = 0 and b = -2; thus, a^2-4*b = 8
	var d gf.GF255e
	r := d.Decode(src)
	zz := d.IsZero()
	d.SqrX(&d, 2)
	d.Lsh(&d, 3)
	d.Add(&d, &gf.GF255e_ONE)
	qr := d.Legendre()

	// If r == 0 then input is not a valid field element encoding
	// and we must return -1 (in such a case, Decode() returned 0,
	// and we have zz == 1 and qr == 0).
	// If r == 1:
	//    If zz == 1, then the point is the neutral. In that case,
	//    we have qr == 0.
	//    If zz == 0, then qr == 1 for valid points, -1 for invalid.
	r = (r - 1) + (-r & (zz - 1) & qr)
	return int(int64(r))
}

// Decode a point from exactly 32 bytes. Returned value is 1 if the
// point could be successfully decoded into a non-neutral group element,
// 0 if it could be successfully decoded as the neutral element, or -1
// if it could not be decoded. If the decoding was not successful, then
// the destination structure is set to the neutral element.
//
// This function is constant-time with regard to the decoded value
// and also with regard to the validity status (timing-based side
// channels do not leak whether the value was found to be a valid point).
//
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral element
//   -1   invalid encoding
func (P *Jq255ePoint) Decode(src []byte) int {
	var e, u, t gf.GF255e

	// Decode the field element as the u coordinate.
	r := u.Decode(src)

	// Compute e^2 = (a^2-4*b)*u^4 - 2*a*u^2 + 1
	// jq255e: a = 0 and b = -2
	e.SqrX(&u, 2)
	e.Lsh(&e, 3)
	e.Add(&e, &gf.GF255e_ONE)
	r &= e.Sqrt(&e)

	// Sqrt() already returns the non-negative square root, so
	// there is no sign adjustment on e and u.

	// t = u^2
	t.Sqr(&u)

	// On failure (r == 0), we replace the coordinates with the
	// neutral element.
	P.e.Select(&e, &gf.GF255e_ONE, r)
	P.z.Set(&gf.GF255e_ONE)
	P.u.Select(&u, &gf.GF255e_ZERO, r)
	P.t.Select(&t, &gf.GF255e_ZERO, r)

	// If r == 0 then the source is invalid and we want to return -1.
	// Otherwise, the point is the neutral if and only if u == 0.
	zz := P.u.IsZero()
	return int(int64((r - 1) + (-r & (1 - zz))))
}

// Test whether a point is the neutral element.
// Returned value is 1 for the neutral, 0 otherwise.
func (P *Jq255ePoint) IsNeutral() int {
	return int(P.u.IsZero())
}

// Test whether this structure (P) represents the same point as the
// provided other structure (Q).
// Returned value is 1 if both points are the same, 0 otherwise.
func (P *Jq255ePoint) Equal(Q *Jq255ePoint) int {
	var t1, t2 gf.GF255e
	t1.Mul(&P.u, &Q.e)
	t2.Mul(&P.e, &Q.u)
	return int(t1.Eq(&t2))
}

// Copy a point structure into another.
// A pointer to this structure is returned.
func (P *Jq255ePoint) Set(Q *Jq255ePoint) *Jq255ePoint {
	P.e.Set(&Q.e)
	P.z.Set(&Q.z)
	P.u.Set(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// If ctl == 1, then copy point Q1 into P.
// If ctl == 0, then copy point Q2 into P.
// ctl MUST be 0 or 1. This is a constant-time selection primitive.
func (P *Jq255ePoint) Select(P1, P2 *Jq255ePoint, ctl int) {
	P.e.Select(&P1.e, &P2.e, uint64(ctl))
	P.z.Select(&P1.z, &P2.z, uint64(ctl))
	P.u.Select(&P1.u, &P2.u, uint64(ctl))
	P.t.Select(&P1.t, &P2.t, uint64(ctl))
}

// Set this point to the sum of the two provided points.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Add(P1, P2 *Jq255ePoint) *Jq255ePoint {
	var n1, n2, n3, n4, n5, n6, n7, n8 gf.GF255e

	// n1 <- E1*E2
	// n2 <- Z1*Z2
	// n3 <- U1*U2
	// n4 <- T1*T2
	n1.Mul(&P1.e, &P2.e)
	n2.Mul(&P1.z, &P2.z)
	n3.Mul(&P1.u, &P2.u)
	n4.Mul(&P1.t, &P2.t)

	// n5 <- (Z1 + T1)*(Z2 + T2) - n2 - n4 = Z1*T2 + Z2*T1
	n5.Add(&P1.z, &P1.t)
	n7.Add(&P2.z, &P2.t)
	n5.Mul(&n5, &n7)
	n5.Sub(&n5, &n2)
	n5.Sub(&n5, &n4)

	// n6 <- (E1 + U1)*(E2 + U2) - n1 - n3 = E1*U2 + E2*U1
	n6.Add(&P1.e, &P1.u)
	n8.Add(&P2.e, &P2.u)
	n6.Mul(&n6, &n8)
	n6.Sub(&n6, &n1)
	n6.Sub(&n6, &n3)

	// n7 <- n2 - (a^2-4*b)*n4 = Z1*Z2 - (a^2-4*b)*T1*T2
	n8.Lsh(&n4, 3)
	n7.Sub(&n2, &n8)

	// E3 <- (n2 + (a^2-4*b)*n4)*(n1 - 2*a*n3) + 2*(a^2-4*b)*n3*n5
	n2.Add(&n2, &n8)
	n3.Mul(&n3, &n5)
	n2.Mul(&n2, &n1)
	n3.Lsh(&n3, 4)
	P.e.Add(&n2, &n3)

	// Z3 <- n7^2
	P.z.Sqr(&n7)

	// T3 <- n6^2
	P.t.Sqr(&n6)

	// U3 <- n6*n7  (could also be computed with ((n6 + n7)^2 - Z3 - T3)/2)
	P.u.Mul(&n6, &n7)

	return P
}

// Set this point to the difference of the two provided points (P1 - P2).
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Sub(P1, P2 *Jq255ePoint) *Jq255ePoint {
	var P2n Jq255ePoint
	P2n.e = P2.e
	P2n.z = P2.z
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.Add(P1, &P2n)
}

// Internal type for a point in extended affine (e, u, t) coordinates.
// This is used to speed up some computations but we do not make it
// public so that the API remains simple.
type jq255ePointAffine struct {
	e, u, t gf.GF255e
}

// Set this point to the sum of the two provided points, the second of
// which being in extended affine (e, u, t) coordinates.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) addMixed(P1 *Jq255ePoint, P2 *jq255ePointAffine) *Jq255ePoint {
	var n1, n2, n3, n4, n5, n6, n7, n8 gf.GF255e

	// n1 <- E1*E2
	// n2 <- Z1*Z2 = Z1
	// n3 <- U1*U2
	// n4 <- T1*T2
	n1.Mul(&P1.e, &P2.e)
	// n2.Mul(&P1.z, &P2.z)
	n3.Mul(&P1.u, &P2.u)
	n4.Mul(&P1.t, &P2.t)

	// n5 <- Z1*T2 + T1
	n5.Mul(&P1.z, &P2.t)
	n5.Add(&n5, &P1.t)

	// n6 <- (E1 + U1)*(E2 + U2) - n1 - n3 = E1*U2 + E2*U1
	n6.Add(&P1.e, &P1.u)
	n8.Add(&P2.e, &P2.u)
	n6.Mul(&n6, &n8)
	n6.Sub(&n6, &n1)
	n6.Sub(&n6, &n3)

	// n7 <- n2 - (a^2-4*b)*n4 = Z1 - (a^2-4*b)*T1*T2
	n8.Lsh(&n4, 3)
	n7.Sub(&P1.z, &n8)

	// E3 <- (n2 + (a^2-4*b)*n4)*(n1 - 2*a*n3) + 2*(a^2-4*b)*n3*n5
	n2.Add(&P1.z, &n8)
	n3.Mul(&n3, &n5)
	n2.Mul(&n2, &n1)
	n3.Lsh(&n3, 4)
	P.e.Add(&n2, &n3)

	// Z3 <- n7^2
	P.z.Sqr(&n7)

	// T3 <- n6^2
	P.t.Sqr(&n6)

	// U3 <- n6*n7  (could also be computed with ((n6 + n7)^2 - Z3 - T3)/2)
	P.u.Mul(&n6, &n7)

	return P
}

// Set this point to the difference of the two provided points, the second of
// which being in extended affine (e, u, t) coordinates.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) subMixed(P1 *Jq255ePoint, P2 *jq255ePointAffine) *Jq255ePoint {
	var P2n jq255ePointAffine
	P2n.e = P2.e
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.addMixed(P1, &P2n)
}

// Set this point (P) to the double of the provided point Q.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Double(Q *Jq255ePoint) *Jq255ePoint {
	var n, x, w, j gf.GF255e

	// Doubling EZUT -> XWJ  (1M+3S)
	n.Sqr(&Q.e)
	w.Sqr(&Q.z)
	x.Sqr(&n)
	w.Lsh(&w, 1)
	j.Mul(&Q.e, &Q.u)
	w.Sub(&w, &n)
	j.Lsh(&j, 1)

	// Conversion XWJ -> EZUT  (3S)
	n.Lsh(&x, 1)
	P.z.Sqr(&w)
	P.t.Sqr(&j)
	P.u.Mul(&j, &w)   // could be ((W + J)^2 - Z - T)/2
	P.e.Sub(&n, &P.z)

	return P
}

// Set this point (P) to (2^n)*Q (i.e. perform n successive doublings).
// This function is constant-time with regard to the point values, but
// not to the number of doublings (n); computation time is proportional
// to n.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) DoubleX(Q *Jq255ePoint, n uint) *Jq255ePoint {
	var k, x, w, j gf.GF255e

	if n == 0 {
		P.Set(Q)
		return P
	}

	// First doubling EZUT -> XWJ  (1M+3S)
	k.Sqr(&Q.e)
	w.Sqr(&Q.z)
	x.Sqr(&k)
	w.Lsh(&w, 1)
	j.Mul(&Q.e, &Q.u)
	w.Sub(&w, &k)
	j.Lsh(&j, 1)

	// n-1 doublings in XWJ coordinates (1M+5S each)
	for n --; n > 0; n -- {
		var k1, k2 gf.GF255e

		k1.Sqr(&w)                     // k1 <- W^2
		k2.Lsh(&x, 1).Sub(&k1, &k2)    // k2 <- k1 - 2*X
		w.Add(&w, &k2).Sqr(&w)         // W' <- (W + k2)^2
		k2.Sqr(&k2)                    // k2 <- k2^2
		w.Sub(&w, &k1).Sub(&w, &k2)    // W' <- W' - k1 - k2
		j.Mul(&j, &w)                  // J <- J*W'
		k1.Sqr(&k1).Lsh(&k1, 1)        // k1 <- 2*k1^2
		w.Sub(&k2, &k1)                // W <- k2 - k1
		x.Sqr(&k2)                     // X <- k2^2
	}

	// Conversion XWJ -> EZUT  (3S)
	k.Lsh(&x, 1)
	P.z.Sqr(&w)
	P.t.Sqr(&j)
	P.u.Mul(&j, &w)   // could be ((W + J)^2 - Z - T)/2
	P.e.Sub(&k, &P.z)

	return P
}

// Set P to the opposite of point Q.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Neg(Q *Jq255ePoint) *Jq255ePoint {
	P.e.Set(&Q.e)
	P.z.Set(&Q.z)
	P.u.Neg(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// jq255e endomorphism:
// ====================
//
// We use one of the cases described by Gallant, Lambert and Vanstone in
// their 2001 article:
//   https://www.iacr.org/archive/crypto2001/21390189.pdf
//
// Modulus is p = 2^255 - 18651. Curve has order 2*r, with:
//   r = 2^254 - 131528281291764213006042413802501683931
//
// Let eta = sqrt(-1) mod p. There are two such roots, we use this one:
// 7656063742463026568679823572395325799027601838558345258426535816504372595438
//
// Let phi(e, u) = (e, eta*u)
// In extended coordinates: phi(e, u, t) = (e, eta*u, -t)  (because eta^2 = -1)
//
// phi() is an endomorphism over our group G or order r. For any group
// element P, phi(P) = mu*P for a constant mu which is a square root of -1
// modulo r. With our choice of eta, we have the following mu:
// 23076176648693837106500022901799924463072024427516564762134831823525232195341
//
// Let k be a scalar (integer modulo r). We decompose k into two half-size
// values k0 and k1 such that k = k0 + k1*mu mod r.
//
// Since r = 1 mod 4, it can be written as the sum of two squares. Let u
// and v such that r = u^2 + v^2; the choice is unique up to permutation and
// change of sign; these values can be found by using Lagrange's algorithm
// on the lattice ((mu, 1), (r, 0)). We choose the following:
//
//   u =  34978546233976132960203755786038370577
//   v = 166506827525740345966246169588540045182
//
// Since (u/v)^2 = -1 mod r, value mu is equal to either u/v or v/u (mod r).
// With our choices, mu = u/v mod r.
//
// It can be verified that:
//   r = u^2 + v^2
//   v + mu*u = 0 mod r
//   -u + mu*v = 0 mod r
//
// Given k, we compute integers c and d as:
//   c = round(k*v / r)
//   d = round(k*u / r)
// Note that c and d are nonnegative. Moreover, c and d fit on 127 bits each.
//
// We then have:
//   k0 = k - d*u - c*v
//   k1 = d*v - c*u
//
// It can be shown (see GLV article) that k0^2 + k1^2 <= u^2 + v^2. Since
// u^2 + v^2 = r < 2^254, this implies that |k0| < 2^127 and |k1| < 2^127.
// Thus, k0 and k1 (which are signed integers) can fit on 128 bits each,
// including their sign bit.
//
//
// Rounded division:
// =================
//
// To compute c and d, we need to do a rounded division. We use the
// fact that r = 2^254 - r0 with r0 < 2^127.
//
// Suppose that we have x = k*u or k*v, and we want to compute y = round(x/r).
// Since r is odd, we have:
//   y = round(x/r) = floor((x + (r-1)/2) / r)
// Let z = x + (r-1)/2. We can split z at the 254-bit index:
//   z = z0 + 2^254*z1
// with 0 <= z0 < 2^254, and z1 >= 0. Since k < r and given the values of u
// and v, the maximum value for z1 is about 2^126.97, i.e it fits on 127 bits.
//
// We thus have:
//   z1*r = z1*2^254 - z1*r0 < z1*2^254 < z
// and
//   (z1+2)*r = z1^254 - z1*r0 + 2*2^254 - 2*r0
//            = (z1+1)*2^254 + (2^254 - (z1+2)*r0)
// Since (z1+2)*r0 < 2^254, we have (z1+2)*r > (z1+1)*2^254 > z.
//
// It follows that the rounded division result is necessarily either z1
// or z1+1. We can thus compute that rounded division with the following
// algorithm:
//
//   Input: integer x such that 0 <= x <= (r-1)*max(u,v)
//   Output: round(x / r)
//     1. z <- x + (r-1)/2
//     2. y <- floor(z / 2^254) + 1
//     3. if y*r > z, then y <- y-1
//     4. return y
//
// The most expensive operation is the product y*r. However, we only need
// the sign of z - y*r. We can do that computation as follows:
//    z - y*r = (y-1)*2^254 + z0 - y*2^254 + y*r0
//            = z0 + y*r0 - 2^254
// We thus need to subtract 1 from y if and only if z0 + y*r0 is strictly
// lower than 2^254. y and r0 both fit on 127 bits each, and z0 is less
// than 2^254; we can thus do that computation over 255 bits.

var jq255eEta = gf.GF255e {
	0xD99E0F1BAA938AEE, 0xA60D864FB30E6336,
	0xE414983FE53688E3, 0x10ED2DB33C69B85F }
var jq255eMinusEta = gf.GF255e {
	0x2661F0E4556C2C37, 0x59F279B04CF19CC9,
	0x1BEB67C01AC9771C, 0x6F12D24CC39647A0 }

// Multiply a point Q by a given scalar n.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Mul(Q *Jq255ePoint, n *Jq255eScalar) *Jq255ePoint {
	// Split input scalar into k0 and k1.
	var k0, k1 [2]uint64
	n.SplitMu(&k0, &k1)

	// Recode k0 and fill M with a copy of Q or -Q, depending on
	// whether k0 is negative or not.
	var sd0 [26]byte
	sg := scalar.Recode5SmallSigned(&sd0, &k0)
	M := *Q
	M.u.CondNeg(&M.u, sg)

	// Initialize the low window with i*M for i in 1..16
	// (point i*Q is stored at offset i-1).
	var win0 [16]Jq255ePoint
	win0[0] = M
	win0[1].Double(&M)
	for i := 3; i <= 15; i += 2 {
		win0[i - 1].Add(&win0[i - 2], &M)
		win0[i].Double(&win0[((i + 1) >> 1) - 1])
	}

	// Recode k1 and remember whether the sign of k1 differs from
	// that of k0.
	var sd1 [26]byte
	sg ^= scalar.Recode5SmallSigned(&sd1, &k1)

	// Apply the endomorphism on all points of the low window to
	// fill the high window. If k0 and k1 have different signs,
	// then an extra negation is performed.
	var endo gf.GF255e
	endo.CondNeg(&jq255eEta, sg)
	var win1 [16]Jq255ePoint
	for i := 0; i < 16; i ++ {
		win1[i].e = win0[i].e
		win1[i].z = win0[i].z
		win1[i].u.Mul(&win0[i].u, &endo)
		win1[i].t.Neg(&win0[i].t)
	}

	// Lookup points corresponding to the top digits, and add them
	// to get the initial value of P. The top digits are both
	// nonnegative.
	P.jq255eLookupWindow(&win0, uint(sd0[25]))
	M.jq255eLookupWindow(&win1, uint(sd1[25]))
	P.Add(P, &M)

	// Process other digits from top to bottom.
	for i := 24; i >= 0; i -- {
		P.DoubleX(P, 5)
		M.jq255eLookupWindow(&win0, uint(sd0[i] & 0x1F))
		M.u.CondNeg(&M.u, uint64(sd0[i] >> 7))
		P.Add(P, &M)
		M.jq255eLookupWindow(&win1, uint(sd1[i] & 0x1F))
		M.u.CondNeg(&M.u, uint64(sd1[i] >> 7))
		P.Add(P, &M)
	}

	return P
}

// Multiply the conventional generator by a given scalar n. This is
// functionally equivalent (but faster) to P.Generator().Mul(&P, n).
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) MulGen(n *Jq255eScalar) *Jq255ePoint {
	// Recode input scalar into 5-bit Booth encoding.
	var sd [52]byte
	n.recode5(&sd)

	// Lookup initial accumulator by using the top digit (which is
	// guaranteed nonnegative).
	var Ma jq255ePointAffine
	Ma.jq255eLookupWindowAffine(&jq255eWin_G195_eut, uint(sd[51]))
	P.e = Ma.e
	P.z = gf.GF255e_ONE
	P.u = Ma.u
	P.t = Ma.t

	// Add points corresponding to top digits of the three other
	// quarter-scalars.
	Ma.jq255eLookupWindowAffine(&jq255eWin_G_eut, uint(sd[12] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[12] >> 7))
	P.addMixed(P, &Ma)
	Ma.jq255eLookupWindowAffine(&jq255eWin_G65_eut, uint(sd[25] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[25] >> 7))
	P.addMixed(P, &Ma)
	Ma.jq255eLookupWindowAffine(&jq255eWin_G130_eut, uint(sd[38] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[38] >> 7))
	P.addMixed(P, &Ma)

	// Process all other digits from high to low. We process the
	// four quarter-scalars in parallel.
	for i := 11; i >= 0; i -- {
		P.DoubleX(P, 5)

		Ma.jq255eLookupWindowAffine(&jq255eWin_G_eut, uint(sd[i] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255eLookupWindowAffine(&jq255eWin_G65_eut, uint(sd[i + 13] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 13] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255eLookupWindowAffine(&jq255eWin_G130_eut, uint(sd[i + 26] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 26] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255eLookupWindowAffine(&jq255eWin_G195_eut, uint(sd[i + 39] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 39] >> 7))
		P.addMixed(P, &Ma)
	}

	return P
}

// Constant-time lookup of a point in a window. Provided window has
// 16 elements. Input offset ('index') is in the 0..16 range. This
// function sets P to a copy of win[index - 1] if index != 0, or to
// the neutral if index == 0.
func (P *Jq255ePoint) jq255eLookupWindow(win *[16]Jq255ePoint, index uint) {
	// Initialize P to all-zeros.
	P.e = gf.GF255e_ZERO
	P.z = gf.GF255e_ZERO
	P.u = gf.GF255e_ZERO
	P.t = gf.GF255e_ZERO

	// Lookup all values.
	for i := 0; i < 16; i ++ {
		m := int64(index) - int64(i + 1)
		mm := ^uint64((m | -m) >> 63)
		P.e.CondOrFrom(&win[i].e, mm)
		P.z.CondOrFrom(&win[i].z, mm)
		P.u.CondOrFrom(&win[i].u, mm)
		P.t.CondOrFrom(&win[i].t, mm)
	}

	// Set E and Z to 1 if the index is zero.
	mz := uint64((int64(index) - 1) >> 63)
	P.e.CondOrFrom(&gf.GF255e_ONE, mz)
	P.z.CondOrFrom(&gf.GF255e_ONE, mz)
}

// Constant-time lookup of a point in a window. This is similar to
// jq255eLookupWindow(), except that this function works on points in
// extended affine (e, u, t) coordinates.
func (P *jq255ePointAffine) jq255eLookupWindowAffine(win *[16]jq255ePointAffine, index uint) {
	// Initialize P to all-zeros (which is the valid representation
	// of the neutral element).
	P.e = gf.GF255e_ZERO
	P.u = gf.GF255e_ZERO
	P.t = gf.GF255e_ZERO

	// Lookup all values.
	for i := 0; i < 16; i ++ {
		m := int64(index) - int64(i + 1)
		mm := ^uint64((m | -m) >> 63)
		P.e.CondOrFrom(&win[i].e, mm)
		P.u.CondOrFrom(&win[i].u, mm)
		P.t.CondOrFrom(&win[i].t, mm)
	}

	// Set E to 1 if the index is zero.
	mz := uint64((int64(index) - 1) >> 63)
	P.e.CondOrFrom(&gf.GF255e_ONE, mz)
}

// Add to point P a point from a window, given an encoded index.
// THIS IS NOT CONSTANT-TIME.
func (P *Jq255ePoint) addFromWindowVartime(win *[16]Jq255ePoint, ej byte) {
	j := int(ej & 0x1F)
	if j != 0 {
		if ej < 0x80 {
			P.Add(P, &win[j - 1])
		} else {
			P.Sub(P, &win[j - 1])
		}
	}
}

// Add to point P a point from a window (affine), given an encoded index.
// THIS IS NOT CONSTANT-TIME.
func (P *Jq255ePoint) addFromWindowAffineVartime(win *[16]jq255ePointAffine, ej byte) {
	j := int(ej & 0x1F)
	if j != 0 {
		if ej < 0x80 {
			P.addMixed(P, &win[j - 1])
		} else {
			P.subMixed(P, &win[j - 1])
		}
	}
}

// Set this structure (P) to Q*c + s*G, for a 128-bit unsigned integer c
// (provided as two 64-bit limbs, in little-endian order), scalar s,
// and conventional generator G.
// A pointer to this structure (P) is returned.
func (P *Jq255ePoint) Mul128AddMulGenVartime(Q *Jq255ePoint, c *[2]uint64, s *Jq255eScalar) *Jq255ePoint {

	// TODO: use w-NAF instead of a fixed 5-bit schedule
	// (on (e,u) coordinates, overhead for a sequence of doublings is
	// low, so w-NAF is likely to be slightly faster than fixed 5-bit
	// windows)

	// Recode c
	var sd0 [26]byte
	scalar.Recode5Small(&sd0, c)

	// Recode s
	var sd1 [52]byte
	scalar.Recode5(&sd1, (*[4]uint64)(s))

	// Initialize the window with i*Q for i in 1..16
	// (point i*Q is stored at offset i-1).
	var win [16]Jq255ePoint
	win[0] = *Q
	win[1].Double(Q)
	for i := 3; i <= 15; i += 2 {
		win[i - 1].Add(&win[i - 2], Q)
		win[i].Double(&win[((i + 1) >> 1) - 1])
	}

	// Initialize P with the top digits (which are all nonnegative).
	j := int(sd0[25])
	if j == 0 {
		*P = jq255eNeutral
	} else {
		*P = win[j - 1]
	}
	P.addFromWindowAffineVartime(&jq255eWin_G_eut, sd1[25])
	P.addFromWindowAffineVartime(&jq255eWin_G130_eut, sd1[51])

	// Process all other digits in top to bottom order.
	for i := 24; i >= 0; i -- {
		P.DoubleX(P, 5)
		P.addFromWindowVartime(&win, sd0[i])
		P.addFromWindowAffineVartime(&jq255eWin_G_eut, sd1[i])
		P.addFromWindowAffineVartime(&jq255eWin_G130_eut, sd1[26 + i])
	}

	return P
}

// Map a sequence of bytes into a curve element. The mapping is not
// injective or surjective, and not uniform among possible outputs;
// however, any given point has only a limited number of possible
// pre-images by the map. A hash-to-curve process can be built on top
// of this map, as follows:
//  - Hash some input data in 64 bytes, with a secure hash function or
//    XOF (e.g. SHAKE).
//  - Split these 64 bytes into two halves, and map each of them to a
//    point with this map.
//  - Add the two points together.
func (P *Jq255ePoint) MapBytes(bb []byte) *Jq255ePoint {
	// Decode source bytes as a field element. This applies modular
	// reduction.
	var f gf.GF255e
	f.DecodeReduce(bb)

	// Set a flag if the value is zero. We use it at the end to
	// set the result to N in that case. This allows us to assume
	// that f != 0 in the rest of the code.
	fz := f.IsZero()

	// Formulas: we map into a point on the dual curve y^2 = x^3 + bb*x,
	// with bb = -4*b:
	//   x1 = f + (1-bb)/(4*f)
	//   x2 = d*(f - (1-bb)/(4*f))
	// with d = sqrt(-1) (same square root as in the endomorphism).
	// Then, at least one of x1, x2 and x1*x2 is a proper x
	// coordinate for a curve point. We use the following:
	//   x1^3 + bb*x1 = (64*f^7 + 16*(3+bb)*f^5
	//                   + 4*(3-2*bb-bb^2)*f^3 + (1-bb)^3*f) / (64*f^4)
	//   x2^3 + bb*x2 = -d*(64*f^7 - 16*(3+bb)*f^5
	//                      + 4*(3-2*bb-bb^2)*f^3 - (1-bb)^3*f) / (64*f^4)
	//   ((x1*x2)^3 + bb*(x1*x2)) = (x1^3 + bb*x1)*(x2^3 + bb*x2)
	// For a properly deterministic result, we use the square root
	// of the numerator whose least significant bit (as an integer
	// in the 0..p-1 range) is zero; this is what Sqrt() returns.

	var f2, f3, f4, f5, f7 gf.GF255e
	f2.Sqr(&f)
	f4.Sqr(&f2)
	f3.Mul(&f, &f2)
	f5.Mul(&f3, &f2)
	f7.Mul(&f5, &f2)

	// Compute x1 and x2, each as a fraction; the denominator is the
	// same for both (4*f). With bb = 8, we have:
	//   x1num  = 4*f^2 - 7
	//   x2num  = d*(4*f^2 + 7)
	//   x12den = 4*f
	var x1num, x2num, x12den gf.GF255e
	x1num.Lsh(&f2, 2)
	x1num.Sub(&x1num, &gf.GF255e_SEVEN)
	x2num.Lsh(&f2, 2)
	x2num.Add(&x2num, &gf.GF255e_SEVEN)
	x2num.Mul(&x2num, &jq255eEta)
	x12den.Lsh(&f, 2)

	// Compute the values of y1^2 and y2^2 (candidates).
	// With bb = 8, numerator polynomials are:
	//   yy1num = 64*f^7 + 176*f^5 - 308*f^3 - 343*f
	//   yy2num = -d*(64*f^7 - 176*f^5 - 308*f^3 + 343*f)
	// for the common denominator:
	//   yyden  = 64*f^4
	// We directly use the square root of the denominator:
	//   y12den = 8*f^2
	// Note that since 1-bb = -7, and 7 and -7 are not quadratic
	// residues in the field, then neither x1 nor x2 can be zero;
	// moreover, x^2 + bb cannot be 0 for any x, since -bb is not a
	// quadratic residue. Therefore, x1^3 + bb*x1 and x2^3 + bb*x2
	// are both non-zero, and the Legendre symbol is then 1 or -1
	// for each.
	var yy1num, yy2num, y12den, tt gf.GF255e
	yy1num.Lsh(&f7, 6)
	yy2num.Set(&yy1num)
	tt.Mul(&f5, &gf.GF255e_HUNDREDSEVENTYSIX)
	yy1num.Add(&yy1num, &tt)
	yy2num.Sub(&yy2num, &tt)
	tt.Mul(&f3, &gf.GF255e_THREEHUNDREDEIGHT)
	yy1num.Sub(&yy1num, &tt)
	yy2num.Sub(&yy2num, &tt)
	tt.Mul(&f, &gf.GF255e_THREEHUNDREDFORTYTHREE)
	yy1num.Sub(&yy1num, &tt)
	yy2num.Add(&yy2num, &tt)
	yy2num.Mul(&yy2num, &jq255eMinusEta)
	y12den.Lsh(&f2, 3)

	// x3 = x1*x2
	// y3^2 = (y1^2)*(y2^2)
	var x3num, x3den, yy3num, y3den gf.GF255e
	x3num.Mul(&x1num, &x2num)
	x3den.Lsh(&f2, 4)
	yy3num.Mul(&yy1num, &yy2num)
	y3den.Lsh(&f4, 6)

	// Get the Legendre symbols for y1^2 and y2^2.
	ls1 := yy1num.Legendre()
	ls2 := yy2num.Legendre()

	// If ls1 == 1, then we use x1 and yy1. Otherwise, if ls1 == -1
	// and ls2 == 1, then we use x2 and yy2. Otherwise, we use x1*x2,
	// and yy1*yy2.
	qr1 := 1 - (ls1 >> 63)
	qr2 := 1 - (ls2 >> 63)

	var xnum, xden, yynum, yden gf.GF255e
	xnum.Select(&x1num, &x2num, qr1)
	xnum.Select(&xnum, &x3num, qr1 | qr2)
	xden.Select(&x12den, &x3den, qr1 | qr2)
	yynum.Select(&yy1num, &yy2num, qr1)
	yynum.Select(&yynum, &yy3num, qr1 | qr2)
	yden.Select(&y12den, &y3den, qr1 | qr2)

	// Extract the square root of yynum; it cannot fail. The Sqrt()
	// function ensures that we get a least significant bit equal to 0.
	var ynum gf.GF255e
	ynum.Sqrt(&yynum)

	// We now have the point in fractional (x,y) coordinates. We
	// compute the coordinate u = x/y.
	var unum, uden gf.GF255e
	unum.Mul(&xnum, &yden)
	uden.Mul(&xden, &ynum)

	// Apply the isogeny theta'_{1/2} to get a point in the proper
	// group on the right curve:
	//   x' = 4*b*u^2
	//   u' = 2*x/(u*(x^2 - bb))
	var xn, xd, un, ud gf.GF255e
	xn.Sqr(&unum).Lsh(&xn, 3).Neg(&xn)
	xd.Sqr(&uden)
	un.Mul(&xnum, &xden).Mul(&un, &uden).Lsh(&un, 1)
	xden.Sqr(&xden).Lsh(&xden, 3)
	ud.Sqr(&xnum).Sub(&ud, &xden).Mul(&ud, &unum)

	// Convert to EZUT coordinates.
	var en, ed gf.GF255e
	xn.Sqr(&xn)
	xd.Sqr(&xd).Lsh(&xd, 1)
	en.Add(&xn, &xd)
	ed.Sub(&xn, &xd)
	var un2, ud2 gf.GF255e
	un2.Sqr(&un)
	ud2.Sqr(&ud)
	P.e.Mul(&en, &ud2)
	P.z.Mul(&ed, &ud2)
	P.u.Mul(&un, &ud).Mul(&P.u, &ed)
	P.t.Mul(&un2, &ed)

	// If the source value f was zero, then all of the above is
	// invalid, and we force the point to the neutral.
	P.e.Select(&gf.GF255e_ONE, &P.e, fz)
	P.z.Select(&gf.GF255e_ONE, &P.z, fz)
	P.u.Select(&gf.GF255e_ZERO, &P.u, fz)
	P.t.Select(&gf.GF255e_ZERO, &P.t, fz)

	return P
}

// =====================================================================
// Precomputed windows in extended affine (e, u, t) coordinates:
//   i*G         for i = 1..16
//   i*2^65*G    for i = 1..16
//   i*2^130*G   for i = 1..16
//   i*2^195*G   for i = 1..16

// Points i*G for i = 1 to 16, extended affine coordinates
var jq255eWin_G_eut = [16]jq255ePointAffine {
	// G * 1
	{
		e: gf.GF255e { 0x0000000000000003, 0x0000000000000000,
		               0x0000000000000000, 0x0000000000000000 },
		u: gf.GF255e { 0x0000000000000001, 0x0000000000000000,
		               0x0000000000000000, 0x0000000000000000 },
		t: gf.GF255e { 0x0000000000000001, 0x0000000000000000,
		               0x0000000000000000, 0x0000000000000000 },
	},
	// G * 2
	{
		e: gf.GF255e { 0xD0FAC687D6342FD1, 0x687D6343EB1A1F58,
		               0x343EB1A1F58D0FAC, 0x1A1F58D0FAC687D6 },
		u: gf.GF255e { 0xB6DB6DB6DB6D97A3, 0xDB6DB6DB6DB6DB6D,
		               0x6DB6DB6DB6DB6DB6, 0x36DB6DB6DB6DB6DB },
		t: gf.GF255e { 0x0A72F05397827791, 0x05397829CBC14E5E,
		               0x829CBC14E5E0A72F, 0x414E5E0A72F05397 },
	},
	// G * 3
	{
		e: gf.GF255e { 0x26AFA803D61A9E2F, 0xAD8273027F0E7D48,
		               0x4D5065C0BA270925, 0x6E6BA44DDB3919FD },
		u: gf.GF255e { 0xC2F21347C4043E79, 0x6B1CEBA6066D4156,
		               0xAB617909A3E20224, 0x12358E75D30336A0 },
		t: gf.GF255e { 0xC4FAF5442BDDB3C7, 0xC58EF652F0485A50,
		               0x0509961D71E284EF, 0x7287BBB2DC59141C },
	},
	// G * 4
	{
		e: gf.GF255e { 0xBCE7BEB9F8390C16, 0xCBFF478EE825BA04,
		               0x96E4BB9C95BDA924, 0x0F1371769561D944 },
		u: gf.GF255e { 0x65A29F71130DB4AD, 0x9F71130DFA47C8BB,
		               0x130DFA47C8BB65A2, 0x7A47C8BB65A29F71 },
		t: gf.GF255e { 0x4D0A213B4402088D, 0x853223D7F44E59F2,
		               0x03ADCBE22101F311, 0x2375E8119918E929 },
	},
	// G * 5
	{
		e: gf.GF255e { 0xD23D2C8BE875C86A, 0x1BD8155773C41197,
		               0x74304444BCDB09C0, 0x3A3E1251980D6493 },
		u: gf.GF255e { 0x1F2B6B08DA5B43EE, 0xE40F8B8BC44A0C63,
		               0x5866F1F8B35FB70C, 0x185034D250F768D7 },
		t: gf.GF255e { 0xC91927493D361051, 0xE00C1E20C1C66FF4,
		               0x8982206A724B43CC, 0x3E3560E7BB5DF4DA },
	},
	// G * 6
	{
		e: gf.GF255e { 0x643AD390229AD5F2, 0xC712074807471E77,
		               0x1673C5B96A6D2E2A, 0x124FA9DF5844C804 },
		u: gf.GF255e { 0x0BD0C5F1F91D6B18, 0xBB4A410D263610A7,
		               0xA1AB0B9D98F35F00, 0x4FA6D8B6AFDDC92B },
		t: gf.GF255e { 0x355D1614AEB11ACD, 0x76ED99CCAEE9D26F,
		               0xD7991971E94A460E, 0x34F3562FDA88753E },
	},
	// G * 7
	{
		e: gf.GF255e { 0xE8D6CBC7678229FB, 0xD86308C5232C88D0,
		               0xC36F69694BFC77E7, 0x4C5C3FA382AAF7AC },
		u: gf.GF255e { 0x7EB52414159EF4EA, 0xB885C9D1EB4CC9E1,
		               0x350914B3EE64BF7F, 0x6DD8CDFA520AED5A },
		t: gf.GF255e { 0x59DAD0E634C75544, 0x818C73930C2A0899,
		               0x0957AB7A60AC1520, 0x56861F4D0A217C1C },
	},
	// G * 8
	{
		e: gf.GF255e { 0x4E11C11353187B28, 0xAD881FFB515E9113,
		               0x86AD97E42FD5BA49, 0x2C922B2CB2A6146C },
		u: gf.GF255e { 0x99DA8C93EB513A8B, 0x0706B8B95DEDFC87,
		               0xC54D8F471F778CE9, 0x4766315BFA2E63E5 },
		t: gf.GF255e { 0x4D9B8F5639729F9A, 0x89B68A9C8B0077A8,
		               0xF3C520B8FCA311FD, 0x532698BCB811270A },
	},
	// G * 9
	{
		e: gf.GF255e { 0x6FB66DF6B52FBCC4, 0x675E5BCC38AA1784,
		               0x55B6D3E8852C1B0B, 0x2289F3ABFA293050 },
		u: gf.GF255e { 0xA84A27A9D0A08E61, 0x27E9084D132CCAC1,
		               0x498C7D8B01F68C40, 0x6957FDFF940E4159 },
		t: gf.GF255e { 0x8D2F2DE6815F2EFF, 0x76CA668F88C812F9,
		               0x56244B8A32B42796, 0x431DA1A672CB2D3C },
	},
	// G * 10
	{
		e: gf.GF255e { 0x8EB44683ACC048BE, 0x803FC1C6AAD5CD46,
		               0xAF5539730762C505, 0x30290CF961B06E3A },
		u: gf.GF255e { 0x3AA366BBB889903E, 0x55838146CC140A37,
		               0x4AA37581A9B6AD5E, 0x7B37113C916F803C },
		t: gf.GF255e { 0xD912EBC4E1C6283E, 0xC70EAC518AE5C163,
		               0x9EDDA370E828C438, 0x252DC97C189ECFD9 },
	},
	// G * 11
	{
		e: gf.GF255e { 0x38A7E99599D93A3B, 0x09DF0EB0A3919A65,
		               0x6F385F29F643AF23, 0x467C84CA2424A548 },
		u: gf.GF255e { 0xA9A8911D864E7F82, 0x65CF6B9CAB741725,
		               0x8C133221E772B327, 0x158521078CD1F209 },
		t: gf.GF255e { 0x41583C9A8F92D685, 0xFAE4DC5553E938EB,
		               0xC3FC1F026C5406EA, 0x5D4A07E9BC1F036B },
	},
	// G * 12
	{
		e: gf.GF255e { 0xA49A13F883BF7B93, 0x06E6D0CF8DC7807E,
		               0x4B3B87C8238DA4B1, 0x66FF84AF5A9F7719 },
		u: gf.GF255e { 0xE4C1C725087640AA, 0xFB902D6A3EF5D5E0,
		               0x53EF35932E1297EB, 0x67E65CF7E1787343 },
		t: gf.GF255e { 0x4203ACE2FF9309B4, 0xAE5BB5318E506208,
		               0x4742F3CB3DEB52CB, 0x2213A3D93959DA85 },
	},
	// G * 13
	{
		e: gf.GF255e { 0xA99622782ED04723, 0x22587604B5D2A716,
		               0x3E7BB13DFFB6AD2D, 0x6E7036BE9D4C885B },
		u: gf.GF255e { 0x90F8839881061965, 0x67D0394FF2BFCB98,
		               0x913200FCCD1396D8, 0x17F96D76306A3580 },
		t: gf.GF255e { 0xF1F0EC984099DC93, 0xE02396E9E43361F5,
		               0x028EBB02AB0AE384, 0x0E2364672DB22F61 },
	},
	// G * 14
	{
		e: gf.GF255e { 0xC6C60BB30F2521F0, 0xA59973DD5CB6D116,
		               0x069708CC706DD30D, 0x74C367ABF8A08989 },
		u: gf.GF255e { 0x05B49673D2AC4172, 0xA016A6890D77E4E6,
		               0x7C6DAA970635E1C0, 0x42C8034547A6A04A },
		t: gf.GF255e { 0x0266FAA875DEF4DC, 0x41B211E505C5A659,
		               0xE13C4A7639E5E234, 0x0C4AC28DE6AF9B7D },
	},
	// G * 15
	{
		e: gf.GF255e { 0x543AA085E86224A7, 0x626226C09EA21055,
		               0x257B5FE5EE7E01D9, 0x1110D92782C497CD },
		u: gf.GF255e { 0x7FFA4AF719120727, 0x705D12571BF74984,
		               0x4AD1FA649FAE1F07, 0x2F4CA2B6265D7456 },
		t: gf.GF255e { 0xB111F1B5F7C6525E, 0x54BD0CFFC1B29AC7,
		               0xCC7CCE327009957D, 0x0CCF7FF00D563132 },
	},
	// G * 16
	{
		e: gf.GF255e { 0x97DA44A024F31B2E, 0xF8FAE043DB5120DD,
		               0x03D9F770D7F5F415, 0x676824C9A296F053 },
		u: gf.GF255e { 0x2D808316E1227049, 0x15064C9132683177,
		               0x706D8A1F41E90ED8, 0x251A19311A6DB76E },
		t: gf.GF255e { 0x95191DCA9E05F91A, 0x0DB49CC10C6EE0A8,
		               0x7C16D8FF7BF95128, 0x2C8D5EC4B15D04AE },
	},
}

// Points i*(2^65)*G for i = 1 to 16, affine extended format
var jq255eWin_G65_eut = [16]jq255ePointAffine {
	// (2^65)*G * 1
	{
		e: gf.GF255e { 0x886FF6278E08211F, 0x4629856EF94734AB,
		               0xAB3BAC1AB41DD08C, 0x6CE0B77B634A53B5 },
		u: gf.GF255e { 0xAD5F8FBFC596FA71, 0x415893549DE223FF,
		               0x395D2181E50A4384, 0x1B313D36A8A7626E },
		t: gf.GF255e { 0x17E71DD1EB9D832B, 0x222B7C0CD599D9CB,
		               0x48E8393CFF13EC0E, 0x6E78594A21A5AAD2 },
	},
	// (2^65)*G * 2
	{
		e: gf.GF255e { 0x111922C72716A209, 0xE4584C1E3417B3E6,
		               0x8AA629F4CE6C2582, 0x7521310F8DD93292 },
		u: gf.GF255e { 0x2B7C2F71889E3F33, 0x4CA4C049A5E65CF4,
		               0x5CD27D909976BFE7, 0x0BE56F359985D602 },
		t: gf.GF255e { 0xC699D8CD1BE2027B, 0xC9233D8219ECC70E,
		               0x6C805A5DABAECD17, 0x40534B306985620E },
	},
	// (2^65)*G * 3
	{
		e: gf.GF255e { 0x33F6746ED6A45A0E, 0x5A40972FCDA92791,
		               0x7B11CD4B6C9CC38D, 0x390940AD8D0885BF },
		u: gf.GF255e { 0x4A504A5DED61CB7F, 0xAF41508342D7801D,
		               0x0519A68AAB4295EB, 0x098D3AB90B09C2B4 },
		t: gf.GF255e { 0xA06F9DAD59A456C6, 0x92D534263FFE78C9,
		               0x0696A39545C181DE, 0x41A590927B5F6FD7 },
	},
	// (2^65)*G * 4
	{
		e: gf.GF255e { 0x6C4E581C32A76153, 0xEADD9BBED6E7A570,
		               0xBD9A0206C0DF4A3B, 0x10D60536FA934AC1 },
		u: gf.GF255e { 0xD98BB32E81B4D89F, 0x5A0FCD3AE1067F08,
		               0x3B845EF3AFAD3191, 0x3540EB32AA6E2C23 },
		t: gf.GF255e { 0x1D5909B1A76336E1, 0x4563C7AD47CE6A39,
		               0x893FFC0CA3CB3C98, 0x6470D11011381CFF },
	},
	// (2^65)*G * 5
	{
		e: gf.GF255e { 0x648443E48591C275, 0x74427A8818227B9E,
		               0x9AE9B0736E87C51E, 0x742236A2B9985165 },
		u: gf.GF255e { 0xDB2EB15DD80EB7AB, 0x695C5863633AC0BB,
		               0x9DD36D925FB45810, 0x4B36B3BE374E74CB },
		t: gf.GF255e { 0x3B02030497F19A1B, 0xBCF6E190CBC99BC3,
		               0x52217A6F3280BAD3, 0x7563AAD0685E3293 },
	},
	// (2^65)*G * 6
	{
		e: gf.GF255e { 0xD8AD585D0C14DC83, 0x748A8C77D6669E2A,
		               0x9318569D6BC1A5D5, 0x637BA56A42C4587D },
		u: gf.GF255e { 0x2473D6DB1425C001, 0x763067186E58108F,
		               0x0C0776C9F6ED32DE, 0x143CF8A090ACB085 },
		t: gf.GF255e { 0x73CDF9E505579C38, 0xDE3B04BA9040AF87,
		               0x3A3345F98418DD26, 0x2C301B6BCC1F09C1 },
	},
	// (2^65)*G * 7
	{
		e: gf.GF255e { 0xE14526AA9D051B56, 0x0E16666D8273B21D,
		               0x1094A46B6F56FA43, 0x5C780EAA96A1BB91 },
		u: gf.GF255e { 0xE95E218127D26209, 0x33088969E8FB612B,
		               0xBAC06821F2BF788C, 0x7EE7B8C1D61C83DD },
		t: gf.GF255e { 0x0474E6971F3502A9, 0x624FFDB35FA59240,
		               0x617E1C89A1A25303, 0x72AB912B16C33383 },
	},
	// (2^65)*G * 8
	{
		e: gf.GF255e { 0xE800A094D67CC3ED, 0xD6C684EDF26076D9,
		               0xE95C514D17DD5343, 0x6B3865BE473732FF },
		u: gf.GF255e { 0x4C033B2301944CB8, 0x3260DA08A5D337CB,
		               0x34DFEFB2A6C39FFC, 0x235872299B5B142E },
		t: gf.GF255e { 0x439324FCD04090FB, 0x08AA04266FC23C15,
		               0x7B1EF0DEBBFAA0EE, 0x2961CB37E34C35A1 },
	},
	// (2^65)*G * 9
	{
		e: gf.GF255e { 0xE7847250A19EF549, 0xDD56EA4894B8F6C1,
		               0xF8815AEC9C0105B4, 0x32138C48916B15C7 },
		u: gf.GF255e { 0x6742B67C6086DF92, 0x53193A718D75331F,
		               0xADD356D1CA14E352, 0x07DF0CC5D7C9E88D },
		t: gf.GF255e { 0x5C9C41C3B8594D9B, 0x643500CDE6FADF0D,
		               0x0D31937E03003B68, 0x713B783B3766078D },
	},
	// (2^65)*G * 10
	{
		e: gf.GF255e { 0x25F90BCEB299A5D0, 0x4FD5981D24E25FBE,
		               0x83A3CE92BEFEC84E, 0x5C9329786A90B198 },
		u: gf.GF255e { 0x7982C54051C4DF08, 0x476A1D01EBBFE2BB,
		               0x3F6064CC7E287D25, 0x79E82E420C17C457 },
		t: gf.GF255e { 0x8CE01EA2E8B54114, 0x73165189E2EF99F8,
		               0x27984BF87B6F0AD1, 0x1813C6ED70DAFA38 },
	},
	// (2^65)*G * 11
	{
		e: gf.GF255e { 0x7CAD33020884167D, 0x383BFB84CDDA791E,
		               0x8D0D667D16879C64, 0x4515BBA1AEAFD937 },
		u: gf.GF255e { 0x00FF128613403E58, 0x854CBABF2D02F041,
		               0xB97069DD65593890, 0x3CC8BC3168B5A376 },
		t: gf.GF255e { 0xC303EBB8F65D98BA, 0xF8F2F98E5BEF4381,
		               0x206DA25D9FE2783E, 0x058DA5FED462C39D },
	},
	// (2^65)*G * 12
	{
		e: gf.GF255e { 0x7118E3577FDDE52D, 0x4F20F40E6D616E5E,
		               0x68592F914D7E5167, 0x3E658213D95C2805 },
		u: gf.GF255e { 0x3688782C6588D284, 0xA50B5985163B51DB,
		               0x23447D42F7311FA5, 0x51CBB2B65B751D79 },
		t: gf.GF255e { 0x877FB541EB5BEC17, 0xDC1EF290D77A0305,
		               0x15EB451E43025B6A, 0x60570F440CB85A8A },
	},
	// (2^65)*G * 13
	{
		e: gf.GF255e { 0xBF9758CD45295FBE, 0x0DAC5E2EBE8D1017,
		               0x3CA545A989E728B3, 0x40255DC81F644B88 },
		u: gf.GF255e { 0x6AF277E67A9FAA68, 0xFDEE555C9009703F,
		               0xAA36C89F8D5B502A, 0x052051EDAFC87131 },
		t: gf.GF255e { 0xDED77C4486AAB04C, 0x9B85D45C9AA3A8E1,
		               0x0F6D233DC4F2C4BA, 0x50C64067319D2237 },
	},
	// (2^65)*G * 14
	{
		e: gf.GF255e { 0x13D6BB104910DE30, 0x52126D4376848E4E,
		               0x8E7CF5FB0130A4FF, 0x631EFE8EF7F8E0F4 },
		u: gf.GF255e { 0x86A494FC08EB56DD, 0xF694DF6ED62219A0,
		               0x214282B022098B86, 0x7182E92FA95620EB },
		t: gf.GF255e { 0x8A96C7C2F4D51B6A, 0x6E06A2EFFBF79614,
		               0x2F8383D036FFDEC3, 0x136A3753275A04FE },
	},
	// (2^65)*G * 15
	{
		e: gf.GF255e { 0x714774712F78D8DC, 0x533027C62A4254F0,
		               0xBF13A51C8C3DB802, 0x09C0518DA671868A },
		u: gf.GF255e { 0xD25E4076E71417FD, 0xAE53B9A6D617E037,
		               0xEDF6131ABB08D620, 0x406E88577206C6C4 },
		t: gf.GF255e { 0xA7B19AB7E5261DBE, 0x5FD6950193780470,
		               0x69CB583FD3DDA8C2, 0x7A35DA9A106BD4A6 },
	},
	// (2^65)*G * 16
	{
		e: gf.GF255e { 0xB7E60EBD3E152D96, 0x16B77AD03CBC9CE8,
		               0xA8D3D5EBA8F300F7, 0x0DA857B86203938C },
		u: gf.GF255e { 0x7D37BBB6D0781C73, 0x1ADF4CFE82DAF6AB,
		               0x6CB90D1DDFBA601C, 0x7B916A823C090D35 },
		t: gf.GF255e { 0xFE1CAB3D1875E47C, 0x75D85076E298321B,
		               0xAA62937A623721E5, 0x78BE6830E382EC23 },
	},
}

// Points i*(2^130)*G for i = 1 to 16, affine extended format
var jq255eWin_G130_eut = [16]jq255ePointAffine {
	// (2^130)*G * 1
	{
		e: gf.GF255e { 0x869B548FDE63606E, 0x04DA93AE9BE27159,
		               0x9FBB8D6BBC2E0657, 0x06CB9DE70E47A525 },
		u: gf.GF255e { 0x6A6F8767CAF58BCB, 0xEF9A2E5ACD520CD9,
		               0x2B998E19EE40437C, 0x1E3A7692F3E02AB1 },
		t: gf.GF255e { 0x12DCF83EC50BD952, 0xC1FEC4C782CB44C6,
		               0x2FCDD9E72158E8D4, 0x6383719F46FAC0AE },
	},
	// (2^130)*G * 2
	{
		e: gf.GF255e { 0xB08E5A163D54F276, 0x137941C5ADF1CAD3,
		               0xE12975DE1B6AA715, 0x69E07759E1492E65 },
		u: gf.GF255e { 0x97643EA7C28259E8, 0x64C33BBEA0416456,
		               0xAC5EBA85AFFFBCEB, 0x1DE0359F8936CEEA },
		t: gf.GF255e { 0xF32ACE061879EA59, 0xC8A5F632AB9427B0,
		               0xF83CBCACE081587F, 0x3A18C76587D69006 },
	},
	// (2^130)*G * 3
	{
		e: gf.GF255e { 0xA69544938C9DC498, 0x9225D9C973152188,
		               0x7A63055985C3FC53, 0x7F9D9B87A431B60A },
		u: gf.GF255e { 0x4A8F6858185A63C2, 0xC29D0227105E6338,
		               0x4FA122A357313D72, 0x62BAF3EF3C842009 },
		t: gf.GF255e { 0x183908457B30E0BC, 0x799295C376EF453F,
		               0xC5FE42DEAF33DE83, 0x54C34BDF628B654A },
	},
	// (2^130)*G * 4
	{
		e: gf.GF255e { 0x607C64DAEC32F903, 0x8777F23BA9BA461C,
		               0x6E7804AC92C3BE8D, 0x180B366A5F951D11 },
		u: gf.GF255e { 0x17625F88FFC35397, 0xEF697901DA783099,
		               0xDCCCD3E1EF458EDC, 0x37EDC360CFBBEDE2 },
		t: gf.GF255e { 0x5F7FC0E7A0B94FB0, 0x3E9CC11638D69625,
		               0x62C64632720B072F, 0x515891835D836913 },
	},
	// (2^130)*G * 5
	{
		e: gf.GF255e { 0xB5DA3AAFCEA3BDDC, 0x75DE8E82AE7AD735,
		               0xA458CD165446291A, 0x552022F75A7FF671 },
		u: gf.GF255e { 0x055D6CCDE9EC95DF, 0xEB86F24ADFCADBF5,
		               0x6DDAC4AB9F7AE17C, 0x282E209409ADD692 },
		t: gf.GF255e { 0x60BA3A7F6B02F394, 0x5AF736BE387EA15D,
		               0x541CE0C53EC38691, 0x775F596632FB79E4 },
	},
	// (2^130)*G * 6
	{
		e: gf.GF255e { 0x7C2346E4D222B565, 0xA45999A2258AEB93,
		               0x66ECF43864D273CD, 0x5EF32CC1B5A5952D },
		u: gf.GF255e { 0xBA9DA77EDCAD7528, 0xB6FA24A8892C93FF,
		               0xF406AD1FFA1ACD6B, 0x7E4F46D5131E195C },
		t: gf.GF255e { 0x85A61A2DEA698452, 0xEB5AD1F6B9C3F150,
		               0xA68B3F0189B1BA23, 0x3306A93BCA2FEDC4 },
	},
	// (2^130)*G * 7
	{
		e: gf.GF255e { 0x70FD1203988BD912, 0xE9BDF6AD5FC76584,
		               0x5E05B4C4DB704F5D, 0x4830813F189FDCFE },
		u: gf.GF255e { 0x9E343B6FEBEF1413, 0xC1293137B7F95D67,
		               0x6FB5672D05CA0B5B, 0x372239482AEC172D },
		t: gf.GF255e { 0xEDC0FD9628A8DD36, 0xF796A60E89BAA4F4,
		               0x21A50BBDCF12747E, 0x35266C9ACED8D90B },
	},
	// (2^130)*G * 8
	{
		e: gf.GF255e { 0x68EE4317D1F7E9CE, 0x27214239C8731082,
		               0x5D86EBC524ED9795, 0x303B08CF4D486C20 },
		u: gf.GF255e { 0x21232B19EA446499, 0x0F110F3BEDDC580B,
		               0x7BFEE5167F928B59, 0x4E1CF6CCC48289F4 },
		t: gf.GF255e { 0x1768F10E04969A8C, 0x0834BAF860B519A9,
		               0x94AAE19ECEC7737D, 0x1DFD019EB5360B34 },
	},
	// (2^130)*G * 9
	{
		e: gf.GF255e { 0xFA3E5E1F454A9789, 0x9259F250CE8DDF8F,
		               0x2830FF4E2FBB9A1B, 0x4262554CF57EA2FC },
		u: gf.GF255e { 0xFFFA0532A28FF940, 0x06134380E791A9D0,
		               0x24120A75934E696A, 0x6F671B022C83BC57 },
		t: gf.GF255e { 0x0BE5ED639180368F, 0xBB03E3FDF6E447C2,
		               0x4CA31E3B83E29BB6, 0x03BBDFC635166005 },
	},
	// (2^130)*G * 10
	{
		e: gf.GF255e { 0xE1B24ECCE93AA158, 0xBC5D29A73202BCE6,
		               0x2970E4C0A4D753B1, 0x776F6C1C18E51DF5 },
		u: gf.GF255e { 0x0EE1B5003B071E37, 0x3EC320448D80BCA7,
		               0x108D5A1C5EC7534A, 0x64ECAAC0EBA1A511 },
		t: gf.GF255e { 0x10E500131535C514, 0xD00F9D216D74A32A,
		               0x8C125FB0A13E6050, 0x4BC50BF7A3AC991A },
	},
	// (2^130)*G * 11
	{
		e: gf.GF255e { 0x10DC6E263EB40DBD, 0x31DC4E881F4F95D8,
		               0x2C33E16E5D45B50A, 0x3622E7BF2A3899E8 },
		u: gf.GF255e { 0x6CC3477AD4AB70A1, 0x0BC92225712BA4D4,
		               0xD51C733D4564D8F4, 0x029CAFA5599B372A },
		t: gf.GF255e { 0x3D63011A6518ADF3, 0x84FD6E87816D4B31,
		               0x396102BA67CD8FA6, 0x30F0D235E547B2AD },
	},
	// (2^130)*G * 12
	{
		e: gf.GF255e { 0x5424663282E719F6, 0xB3E38AE4779DEFF8,
		               0x3F2A874D3EE0C99A, 0x70F7A39B8561373F },
		u: gf.GF255e { 0xCB9E1D5CB7854025, 0x960067A5C064493E,
		               0x41F420B5A91ED717, 0x6E81F9FA7E661D8D },
		t: gf.GF255e { 0x96737B3BD1B3FF89, 0xFA5B99758B260B44,
		               0x4E40C5DCD5EB0507, 0x797141D1E37DDC45 },
	},
	// (2^130)*G * 13
	{
		e: gf.GF255e { 0x3A2AA626DE84F742, 0xC52014D6E24C4D0F,
		               0x87620DD6665947AC, 0x4489AC6538239FA8 },
		u: gf.GF255e { 0x1A722773E489DDB8, 0xF94983893D4AABD6,
		               0x59F3D4C5BB3DFDCC, 0x653EE371D2801E6A },
		t: gf.GF255e { 0xA14B344A108032A6, 0x336E96DD99975786,
		               0x3AF72BF16ED6198C, 0x6E8DE13723DFA5BC },
	},
	// (2^130)*G * 14
	{
		e: gf.GF255e { 0x5BDDDC0121EA7956, 0x30C254D3EE996605,
		               0xE23DA12DBE2A729B, 0x0698AEC8D7B177EE },
		u: gf.GF255e { 0x460E164D09693F50, 0x96FABF7744D22EC2,
		               0x216A1928595E868E, 0x50E1BEE9AC402680 },
		t: gf.GF255e { 0x2EA3B4425FE17CBC, 0x3076D3BE8227BB81,
		               0x73999AF999779B03, 0x52BC8B51287FBDD0 },
	},
	// (2^130)*G * 15
	{
		e: gf.GF255e { 0x4DE169D233632097, 0x7269B6542711C087,
		               0xC7562833544AC635, 0x67C8ED8C09C7A4DA },
		u: gf.GF255e { 0x442C914C64C6EE61, 0x5486463AAA3D41CD,
		               0x2323BA05744FB271, 0x5CE94782B63D2983 },
		t: gf.GF255e { 0x5E85E0F841CFEA05, 0xFE575987C8449D15,
		               0x4B8F046B40C3632A, 0x79B75334C85A090C },
	},
	// (2^130)*G * 16
	{
		e: gf.GF255e { 0xA1638BEC45B50B50, 0xB956B5A6669E52E3,
		               0xAFF58E0E6F53165A, 0x5F00BEB6EDB8A088 },
		u: gf.GF255e { 0x20DDE2D9560BD063, 0x68337F979386B815,
		               0x9CAE33A6B5F9B94C, 0x0F2ED8418B17674E },
		t: gf.GF255e { 0x42082E618690FF50, 0x3721E53E5901899E,
		               0xBB88653D342DE052, 0x2EED8F30CF10FA1C },
	},
}

// Points i*(2^195)*G for i = 1 to 16, affine extended format
var jq255eWin_G195_eut = [16]jq255ePointAffine {
	// (2^195)*G * 1
	{
		e: gf.GF255e { 0x259381D52FD3D48B, 0x05A13379CA5B1805,
		               0x4F9CCA6483A60AC4, 0x26FE2B2D2FF7C515 },
		u: gf.GF255e { 0xC6C4EA52864B8022, 0x3FACF03027F2FE05,
		               0x5A78F8FDAFE0F2B2, 0x7A2058682117A352 },
		t: gf.GF255e { 0x827479FBF869915D, 0xD2369B352A0FAC70,
		               0x759EDBD5AA299C4A, 0x6CF49CB73C1C85C0 },
	},
	// (2^195)*G * 2
	{
		e: gf.GF255e { 0xBEA01D24F6BA6D06, 0x2F4C099F8F9811D4,
		               0xFEFB94253EE54433, 0x26E11FA33F986A71 },
		u: gf.GF255e { 0xE7A9C455FDCCE69F, 0xB043C24E23C52866,
		               0xCBC1DD8A3179B032, 0x597FE7EC4E366C38 },
		t: gf.GF255e { 0xF41D8A2928BB5A33, 0xB52F9C48D79CEDB6,
		               0x31EC395A62DF9B38, 0x5D82A3622AECDB76 },
	},
	// (2^195)*G * 3
	{
		e: gf.GF255e { 0x9902C1E090F2D551, 0xE3C20E9FA8811C74,
		               0x5D9FD304B83A88A9, 0x7A97A73E152AC3D0 },
		u: gf.GF255e { 0xEBF1C192782AD7E7, 0xDAC867CE0228990B,
		               0xB0C0AEB839C2A9BB, 0x5D529C2B2E3222F2 },
		t: gf.GF255e { 0xDF8E63928F1AE0C2, 0x8692B8050BB45ACF,
		               0x4CC66CF8E7017825, 0x25F396E870747870 },
	},
	// (2^195)*G * 4
	{
		e: gf.GF255e { 0x825E9E80D9017500, 0x532E7F73C34493C0,
		               0x23A615DE55EB285B, 0x68BEB6AD6278C4C5 },
		u: gf.GF255e { 0xE5AD1FDA050CE7CF, 0xD5179DCB398FE9E2,
		               0x880F0F9CA2B23DE9, 0x73E9DA1D7C583AB6 },
		t: gf.GF255e { 0x1BED8C4A161AD03A, 0x56A385AC631A3736,
		               0x55CA5A73E2FDED3F, 0x2F2A10845751514B },
	},
	// (2^195)*G * 5
	{
		e: gf.GF255e { 0xEF721415A10674B8, 0x39A98FB9520B23D9,
		               0xDD94B1823583A50F, 0x1A980F7E359D5D64 },
		u: gf.GF255e { 0x81533FD0CED00FE0, 0x2B41B323457375B0,
		               0x3428954D0B0B6412, 0x3FB05C6B656FCDE7 },
		t: gf.GF255e { 0xD6A2CFBECD1FF35F, 0x3EB933A63E59FA2B,
		               0x0156D1B1D6CF146F, 0x1FA8E20753FBE8B0 },
	},
	// (2^195)*G * 6
	{
		e: gf.GF255e { 0xC32013800A7FA85F, 0xCE6BBA92A647E0CA,
		               0xD1328954C1969616, 0x53505DC8D146C39F },
		u: gf.GF255e { 0xABECE8DEAA4DEFF3, 0xA6B25F5370FF8BED,
		               0xC70C1F018B95875D, 0x1EEE7F4380019FB8 },
		t: gf.GF255e { 0xCC3FC741CC1E0562, 0x93B664F242B13EF8,
		               0x5D617DE816798EB4, 0x35C68ADE3CAC8CA8 },
	},
	// (2^195)*G * 7
	{
		e: gf.GF255e { 0x293F3FD57745E14A, 0x023E52D8EE207D6F,
		               0xC95A3B2FA7918F82, 0x203E41FE3594BCCC },
		u: gf.GF255e { 0xC2513D53CE5A6CD2, 0x8AF5B5BD4C9ADB58,
		               0xDF748C1856292D78, 0x1C54D437C147EB47 },
		t: gf.GF255e { 0x02C3EA61A5ABCB56, 0xB56CC897BA7BA956,
		               0xB8E346F880AA5525, 0x791ECB5E7F925E67 },
	},
	// (2^195)*G * 8
	{
		e: gf.GF255e { 0x9E436CB767A86CB5, 0x9EAB45BEF1F6D3FE,
		               0xA2654498E4A81FA1, 0x5A50303BB7FBF10D },
		u: gf.GF255e { 0x961441BFA4853698, 0x76396E28425429D3,
		               0x23187D9C49399AB4, 0x47EC89341C754A72 },
		t: gf.GF255e { 0xB5B36F776AF04917, 0x4E7E59858CDEC91B,
		               0x47DBD9D756A70427, 0x00D090A5B2E0E163 },
	},
	// (2^195)*G * 9
	{
		e: gf.GF255e { 0x1390F65E3B6C3E84, 0xD1492FB1E8EDE015,
		               0x3ADC11E0F4D52DD6, 0x4E1EFB691EF8D7C0 },
		u: gf.GF255e { 0x790DDCB656DA5EBB, 0x899CA30F8A6C1157,
		               0xB055E943E160FF52, 0x0C4E4B67E97A3F02 },
		t: gf.GF255e { 0xE665E42E421C783F, 0x90F6162E6D8BF1FD,
		               0x0E7DEA665667BD29, 0x5C4551BBBCA04267 },
	},
	// (2^195)*G * 10
	{
		e: gf.GF255e { 0x28CA2A1B76E757EA, 0x82E93F7DB218D2B4,
		               0x433AC8B15317BE46, 0x0BBAE0E5BEFBB5CF },
		u: gf.GF255e { 0x7CBAD6A204D1F6B9, 0xEB725D50999A4399,
		               0x7C80807D104B0670, 0x44B8942C5DF07889 },
		t: gf.GF255e { 0x617C6A396D315D12, 0x4555C3E786404CA2,
		               0x4584662C9819C28B, 0x3B052B97144BDDCC },
	},
	// (2^195)*G * 11
	{
		e: gf.GF255e { 0x0D101C917FFDE613, 0x200EE9E2631EC7C4,
		               0xD28F560EA8A29B73, 0x2CCB2CF7E40BF4D0 },
		u: gf.GF255e { 0x30FE96716E7DF796, 0xB6A214969A98317F,
		               0xEB5423DB7543C3F6, 0x7DD81F0AD475BB65 },
		t: gf.GF255e { 0x6E9D90B85FC133E0, 0x4ED760314652F82F,
		               0x76EEE28489A2673B, 0x69561B43851B3032 },
	},
	// (2^195)*G * 12
	{
		e: gf.GF255e { 0x50C89E4AB1499911, 0xA66737D8513DFDB2,
		               0xB3A61F723EE958E7, 0x2DFA0BA127D42687 },
		u: gf.GF255e { 0x041881DC4BB15593, 0x0620628690F070A8,
		               0xF6647FFFA6239BFD, 0x60406418E0A8E484 },
		t: gf.GF255e { 0xB186E710CF384958, 0xCC5BCAEE53A27045,
		               0x5F7FF4B10823ABA8, 0x3C01E53D24508710 },
	},
	// (2^195)*G * 13
	{
		e: gf.GF255e { 0xA54981754921B43D, 0x039A3283169BCBD8,
		               0x16210166F7E6E8DA, 0x0B0051552C1C58B5 },
		u: gf.GF255e { 0x02A5C58CE250DF40, 0x80A9A62960313C1E,
		               0xEDC81102A9A286D9, 0x03C8B0610E5DE932 },
		t: gf.GF255e { 0x9B656127329A1F5C, 0x7F61E2E190853286,
		               0x220189A7901D370E, 0x5C2FDF1BE72A0992 },
	},
	// (2^195)*G * 14
	{
		e: gf.GF255e { 0x1725A4D445249491, 0x117438E7D83A2644,
		               0xB7E063EA06740CDF, 0x6021348E9D2BB8D8 },
		u: gf.GF255e { 0x2359C265188EE74D, 0x3498D65BDB7611FC,
		               0x97D9FFD6286A6BD1, 0x63F775224E36165F },
		t: gf.GF255e { 0xD5B3AA1D5B583047, 0xA31076307FE06CED,
		               0x6492A9AFA98D71BD, 0x534AFBC1B782F441 },
	},
	// (2^195)*G * 15
	{
		e: gf.GF255e { 0x66C6224B1937FD6E, 0x07778715A12C63DB,
		               0xD35D4667B3E6947A, 0x6C435D4A27D08486 },
		u: gf.GF255e { 0x14BD04EFD34FC573, 0x58F89E267B42BEA3,
		               0xDD4D47FA083CC9BF, 0x5C69FCC38DA29629 },
		t: gf.GF255e { 0x98A6F2B6623C4605, 0x4F944F4CDD9551F4,
		               0x5A90BF07D1C9E81B, 0x0B9C09556A326820 },
	},
	// (2^195)*G * 16
	{
		e: gf.GF255e { 0xF4448E39173AD989, 0xB4BD88C5AC6C7FE8,
		               0x6217606AC5C899F0, 0x4DA95AF12CBA888E },
		u: gf.GF255e { 0x805E028481B3D10D, 0x3E6A069F39FCEFCB,
		               0xDA636B907FED771B, 0x162581D9B675A4E1 },
		t: gf.GF255e { 0x52668A8902F702CE, 0x44FDFE7C39061B8D,
		               0xEEAA462F252FD554, 0x02FA9D8FB083E563 },
	},
}
