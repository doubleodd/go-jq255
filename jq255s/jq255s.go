package jq255s

import (
	gf "github.com/doubleodd/go-jq255/internal/field"
	"github.com/doubleodd/go-jq255/internal/scalar"
)

// This file implements operations on curve points for jq255s
// (specifically, on elements of the prime order group defined over
// jq255s). It also provides function for computing with scalar
// values (integers modulo the group order).
//
// API: a point is represented in memory by a Jq255sPoint structure.
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
// A scalar is similarly represented in memory by a Jq255sScalar
// mutable structure. Scalars can be encoded to, and decoded from,
// 32 bytes. Encoding is canonical: a given scalar value can be encoded
// in a unique way, and only one sequence of bytes may decode to that
// scalar.
//
// Unless explcitly documented, all functions here are constant-time.

// Jq255sPoint is the type for a jq255s point.
//
// Default value for a point structure is not valid. The NewJq255sPoint()
// function makes sure to return only initialized structures. If allocating
// a point structure manually, make sure to properly set it to a valid point
// before using it as source.
type Jq255sPoint struct {
	// Internally, we use extended (e,u) coordinates, which have
	// complete and efficient formulas.
	e, z, u, t gf.GF255s
}

// Preallocated neutral point. Do not modify.
var jq255sNeutral = Jq255sPoint {
	e: gf.GF255s { 1, 0, 0, 0 },
	z: gf.GF255s { 1, 0, 0, 0 },
	u: gf.GF255s { 0, 0, 0, 0 },
	t: gf.GF255s { 0, 0, 0, 0 },
}

// Preallocated conventional generator point. Do not modify.
var jq255sGenerator = Jq255sPoint {
	e: gf.GF255s { 0x104220CDA2789410, 0x6D7386B2348CC437,
		       0x55E452A64612D10E, 0x0F520B1BA747ADAC },
	z: gf.GF255s { 1, 0, 0, 0 },
	u: gf.GF255s { 3, 0, 0, 0 },
	t: gf.GF255s { 9, 0, 0, 0 },
}

// Create a new point. The point is set to the group neutral element (N).
func NewJq255sPoint() *Jq255sPoint {
	P := new(Jq255sPoint)
	*P = jq255sNeutral
	return P
}

// Set the point P to the neutral element (N).
// A pointer to this structure is returned.
func (P *Jq255sPoint) Neutral() *Jq255sPoint {
	*P = jq255sNeutral
	return P
}

// Set the point P to the conventional generator (G).
// A pointer to this structure is returned.
func (P *Jq255sPoint) Generator() *Jq255sPoint {
	*P = jq255sGenerator
	return P
}

// Encode a point into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (P *Jq255sPoint) Encode(dst []byte) []byte {
        // Encoded value is u or -u: if e is "negative" (least significant
        // bit of the encoding of e is 1), then we use -u, otherwise
        // we use u.
        var iz gf.GF255s
        iz.Inv(&P.z)
        var r gf.GF255s
        r.Mul(&iz, &P.e)
        nn := r.IsNegative()
        r.Mul(&iz, &P.u)
        r.CondNeg(&r, nn)
        return r.Encode(dst)
}

// Encode a point into exactly 32 bytes.
func (P *Jq255sPoint) Bytes() [32]byte {
	var d [32]byte
	P.Encode(d[:0])
	return d
}

// Test whether a given chunk of 32 bytes is a valid representation of
// jq255s group element. This is faster than actually decoding it.
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral point N
//   -1   invalid encoding
func Jq255sCheckPoint(src []byte) int {
	// The bytes can be decoded as a valid point if and only if all of
	// the following hold:
	//  - The bytes can be decoded as a field element u.
	//  - The value (a^2-4*b)*u^4 - 2*a*u^2 + 1 is a square.
	// On jq255s, a = -1 and b = 1/2; thus, a^2-4*b = -1
	var c, d gf.GF255s
	r := c.Decode(src)
	zz := c.IsZero()
	c.Sqr(&c)
	d.Sqr(&c)
	c.Lsh(&c, 1).Sub(&c, &d).Add(&c, &gf.GF255s_ONE)
	qr := c.Legendre()

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
// 0 if it could be successfully decoded as the neutral point N, or -1
// if it could not be decoded. If the decoding was not successful, then
// the destination structure is set to the neutral N.
//
// This function is constant-time with regard to the decoded value
// and also with regard to the validity status (timing-based side
// channels do not leak whether the value was found to be a valid point).
//
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral point N
//   -1   invalid encoding
func (P *Jq255sPoint) Decode(src []byte) int {
	var e, u, t gf.GF255s

	// Decode the field element as the u coordinate.
	r := u.Decode(src)

	// t = u^2
	t.Sqr(&u)

	// Compute e^2 = (a^2-4*b)*u^4 - 2*a*u^2 + 1
	// jq255s: a = -1 and b = 1/2
	var k gf.GF255s
	k.Lsh(&t, 1)
	e.Sqr(&t).Sub(&k, &e).Add(&e, &gf.GF255s_ONE)
	r &= e.Sqrt(&e)

	// Sqrt() already returns the non-negative square root, so
	// there is no sign adjustment on e and u.

	// On failure (r == 0), we replace the coordinates with the
	// neutral element.
	P.e.Select(&e, &gf.GF255s_ONE, r)
	P.z.Set(&gf.GF255s_ONE)
	P.u.Select(&u, &gf.GF255s_ZERO, r)
	P.t.Select(&t, &gf.GF255s_ZERO, r)

	// If r == 0 then the source is invalud and we want to return -1.
	// Otherwise, the point is the neutral if and only if u == 0.
	zz := P.u.IsZero()
	return int(int64((r - 1) + (-r & (1 - zz))))
}

// Test whether a point is the neutral element N.
// Returned value is 1 for the neutral, 0 otherwise.
func (P *Jq255sPoint) IsNeutral() int {
	return int(P.u.IsZero())
}

// Test whether this structure (P) represents the same point as the
// provided other structure (Q).
// Returned value is true if both points are the same, false otherwise.
func (P *Jq255sPoint) Equal(Q *Jq255sPoint) int {
	var t1, t2 gf.GF255s
	t1.Mul(&P.u, &Q.e)
	t2.Mul(&P.e, &Q.u)
	return int(t1.Eq(&t2))
}

// Copy a point structure into another.
// A pointer to this structure is returned.
func (P *Jq255sPoint) Set(Q *Jq255sPoint) *Jq255sPoint {
	P.e.Set(&Q.e)
	P.z.Set(&Q.z)
	P.u.Set(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// If ctl == 1, then copy point Q1 into P.
// If ctl == 0, then copy point Q2 into P.
// ctl MUST be 0 or 1. This is a constant-time selection primitive.
func (P *Jq255sPoint) Select(P1, P2 *Jq255sPoint, ctl int) {
	P.e.Select(&P1.e, &P2.e, uint64(ctl))
	P.z.Select(&P1.z, &P2.z, uint64(ctl))
	P.u.Select(&P1.u, &P2.u, uint64(ctl))
	P.t.Select(&P1.t, &P2.t, uint64(ctl))
}

// Set this point to the sum of the two provided points.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) Add(P1, P2 *Jq255sPoint) *Jq255sPoint {
	var n1, n2, n3, n4, n5, n6, n7 gf.GF255s

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
	n7.Add(&P2.e, &P2.u)
	n6.Mul(&n6, &n7)
	n6.Sub(&n6, &n1)
	n6.Sub(&n6, &n3)

	// n7 <- n2 - (a^2-4*b)*n4 = Z1*Z2 - (a^2-4*b)*T1*T2
	n7.Add(&n2, &n4)

	// E3 <- (n2 + (a^2-4*b)*n4)*(n1 - 2*a*n3) + 2*(a^2-4*b)*n3*n5
	n3.Lsh(&n3, 1)
	n2.Sub(&n2, &n4)
	n1.Add(&n1, &n3)
	n5.Mul(&n3, &n5)
	n2.Mul(&n2, &n1)
	P.e.Sub(&n2, &n5)

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
func (P *Jq255sPoint) Sub(P1, P2 *Jq255sPoint) *Jq255sPoint {
	var P2n Jq255sPoint
	P2n.e = P2.e
	P2n.z = P2.z
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.Add(P1, &P2n)
}

// Internal type for a point in extended affine (e, u, t) coordinates.
// This is used to speed up some computations but we do not make it
// public so that the API remains simple.
type jq255sPointAffine struct {
	e, u, t gf.GF255s
}

// Set this point to the sum of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) addMixed(P1 *Jq255sPoint, P2 *jq255sPointAffine) *Jq255sPoint {
	var n1, n2, n3, n4, n5, n6, n7 gf.GF255s

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
	n7.Add(&P2.e, &P2.u)
	n6.Mul(&n6, &n7)
	n6.Sub(&n6, &n1)
	n6.Sub(&n6, &n3)

	// n7 <- n2 - (a^2-4*b)*n4 = Z1 - (a^2-4*b)*T1*T2
	n7.Add(&P1.z, &n4)

	// E3 <- (n2 + (a^2-4*b)*n4)*(n1 - 2*a*n3) + 2*(a^2-4*b)*n3*n5
	n3.Lsh(&n3, 1)
	n2.Sub(&P1.z, &n4)
	n1.Add(&n1, &n3)
	n5.Mul(&n3, &n5)
	n2.Mul(&n2, &n1)
	P.e.Sub(&n2, &n5)

	// Z3 <- n7^2
	P.z.Sqr(&n7)

	// T3 <- n6^2
	P.t.Sqr(&n6)

	// U3 <- n6*n7  (could also be computed with ((n6 + n7)^2 - Z3 - T3)/2)
	P.u.Mul(&n6, &n7)

	return P
}

// Set this point to the difference of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) subMixed(P1 *Jq255sPoint, P2 *jq255sPointAffine) *Jq255sPoint {
	var P2n jq255sPointAffine
	P2n.e = P2.e
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.addMixed(P1, &P2n)
}

// Set this point (P) to the double of the provided point Q.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) Double(Q *Jq255sPoint) *Jq255sPoint {
	var k, x, w, j gf.GF255s

	// Doubling EZUT -> XWJ  (1M+3S)
	k.Sqr(&Q.u)
	x.Sqr(&k).Lsh(&x, 3)
	k.Lsh(&k, 1)
	w.Add(&Q.t, &Q.z).Sqr(&w)
	j.Mul(&Q.e, &Q.u).Lsh(&j, 1)
	w.Sub(&k, &w)

	// Conversion XWJ -> EZUT  (3S)
	x.Lsh(&x, 1)
	P.z.Sqr(&w)
	P.t.Sqr(&j)
	x.Sub(&x, &P.z)
	P.u.Mul(&j, &w)   // could be ((W + J)^2 - Z - T)/2
	P.e.Sub(&x, &P.t)

	return P
}

// Set this point (P) to (2^n)*Q (i.e. perform n successive doublings).
// This function is constant-time with regard to the point values, but
// not to the number of doublings (n); computation time is proportional
// to n.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) DoubleX(Q *Jq255sPoint, n uint) *Jq255sPoint {
	var k, x, w, j gf.GF255s

	if n == 0 {
		P.Set(Q)
		return P
	}

	// First doubling EZUT -> XWJ  (1M+3S)
	k.Sqr(&Q.u)
	x.Sqr(&k).Lsh(&x, 3)
	k.Lsh(&k, 1)
	w.Add(&Q.t, &Q.z).Sqr(&w)
	j.Mul(&Q.e, &Q.u).Lsh(&j, 1)
	w.Sub(&k, &w)

	// n-1 doublings in XWJ coordinates (2M+4S each)
	for n --; n > 0; n -- {
		var k1, k2, k3 gf.GF255s

		k1.Mul(&w, &j)                   // k1 <- W*J
		k2.Sqr(&k1)                      // k2 <- k1^2
		k1.Lsh(&k1, 1)
		k3.Add(&w, &j).Sqr(&k3)
		k3.Sub(&k3, &k1)                 // k3 <- W^2 + J^2
		j.Lsh(&x, 1)
		x.Sqr(&k2).Lsh(&x, 3)            // X' <- 8*k2^2
		j.Sub(&j, &k3)
		k3.Sqr(&k3)
		w.Lsh(&k2, 1)
		j.Mul(&j, &k1)                   // J' <- 2*k1*(2*X - k3)
		w.Sub(&w, &k3)                   // W' <- 2*k2 - k3^2
	}

	// Conversion XWJ -> EZUT  (3S)
	x.Lsh(&x, 1)
	P.z.Sqr(&w)
	P.t.Sqr(&j)
	x.Sub(&x, &P.z)
	P.u.Mul(&j, &w)   // could be ((W + J)^2 - Z - T)/2
	P.e.Sub(&x, &P.t)

	return P
}

// Set P to the opposite of point Q.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) Neg(Q *Jq255sPoint) *Jq255sPoint {
	P.e.Set(&Q.e)
	P.z.Set(&Q.z)
	P.u.Neg(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// Multiply a point Q by a given scalar n.
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) Mul(Q *Jq255sPoint, n *Jq255sScalar) *Jq255sPoint {
	// Recode input scalar in 52 signed digits in the -15..+16 range.
	var sd [52]byte
	n.recode5(&sd)

	// Initialize a window with i*Q for i in 1..16
	var win [16]Jq255sPoint
	win[0] = *Q
	win[1].Double(Q)
	for i := 3; i <= 15; i += 2 {
		win[i - 1].Add(&win[i - 2], Q)
		win[i].Double(&win[((i + 1) >> 1) - 1])
	}

	// Use the top digit to initialize the accumulator point M.
	P.jq255sLookupWindow(&win, uint(sd[51]))

	// Process other digits from top to bottom.
	for i := 50; i >= 0; i -- {
		P.DoubleX(P, 5)
		var M Jq255sPoint
		M.jq255sLookupWindow(&win, uint(sd[i] & 0x1F))
		M.u.CondNeg(&M.u, uint64(sd[i] >> 7))
		P.Add(P, &M)
	}

	return P
}

// Multiply the conventional generator by a given scalar n. This is
// functionally equivalent (but faster) to P.Generator().Mul(&P, n).
// A pointer to this structure (P) is returned.
func (P *Jq255sPoint) MulGen(n *Jq255sScalar) *Jq255sPoint {
	// Recode input scalar into 5-bit Booth encoding.
	var sd [52]byte
	n.recode5(&sd)

	// Lookup initial accumulator by using the top digit (which is
	// guaranteed nonnegative).
	var Ma jq255sPointAffine
	Ma.jq255sLookupWindowAffine(&jq255sWin_G195_eut, uint(sd[51]))
	P.e = Ma.e
	P.z = gf.GF255s_ONE
	P.u = Ma.u
	P.t = Ma.t

	// Add points corresponding to top digits of the three other
	// quarter-scalars.
	Ma.jq255sLookupWindowAffine(&jq255sWin_G_eut, uint(sd[12] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[12] >> 7))
	P.addMixed(P, &Ma)
	Ma.jq255sLookupWindowAffine(&jq255sWin_G65_eut, uint(sd[25] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[25] >> 7))
	P.addMixed(P, &Ma)
	Ma.jq255sLookupWindowAffine(&jq255sWin_G130_eut, uint(sd[38] & 0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[38] >> 7))
	P.addMixed(P, &Ma)

	// Process all other digits from high to low. We process the
	// four quarter-scalars in parallel.
	for i := 11; i >= 0; i -- {
		P.DoubleX(P, 5)

		Ma.jq255sLookupWindowAffine(&jq255sWin_G_eut, uint(sd[i] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255sLookupWindowAffine(&jq255sWin_G65_eut, uint(sd[i + 13] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 13] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255sLookupWindowAffine(&jq255sWin_G130_eut, uint(sd[i + 26] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 26] >> 7))
		P.addMixed(P, &Ma)
		Ma.jq255sLookupWindowAffine(&jq255sWin_G195_eut, uint(sd[i + 39] & 0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i + 39] >> 7))
		P.addMixed(P, &Ma)
	}

	return P
}

// Constant-time lookup of a point in a window. Provided window has
// 16 elements. Input offset ('index') is in the 0..16 range. This
// function sets P to a copy of win[index - 1] if index != 0, or to
// the neutral if index == 0.
func (P *Jq255sPoint) jq255sLookupWindow(win *[16]Jq255sPoint, index uint) {
	// Initialize P to all-zeros.
	P.e = gf.GF255s_ZERO
	P.z = gf.GF255s_ZERO
	P.u = gf.GF255s_ZERO
	P.t = gf.GF255s_ZERO

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
	P.e.CondOrFrom(&gf.GF255s_ONE, mz)
	P.z.CondOrFrom(&gf.GF255s_ONE, mz)
}

// Constant-time lookup of a point in a window. This is similar to
// jq255sLookupWindow(), except that this function works on points in
// affine (x, u) coordinates.
func (P *jq255sPointAffine) jq255sLookupWindowAffine(win *[16]jq255sPointAffine, index uint) {
	// Initialize P to all-zeros (which is the valid representation
	// of the neutral element).
	P.e = gf.GF255s_ZERO
	P.u = gf.GF255s_ZERO
	P.t = gf.GF255s_ZERO

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
	P.e.CondOrFrom(&gf.GF255s_ONE, mz)
}

// Add to point P a point from a window, given an encoded index.
// THIS IS NOT CONSTANT-TIME.
func (P *Jq255sPoint) addFromWindowVartime(win *[16]Jq255sPoint, ej byte) {
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
func (P *Jq255sPoint) addFromWindowAffineVartime(win *[16]jq255sPointAffine, ej byte) {
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
func (P *Jq255sPoint) Mul128AddMulGenVartime(Q *Jq255sPoint, c *[2]uint64, s *Jq255sScalar) *Jq255sPoint {

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
	var win [16]Jq255sPoint
	win[0] = *Q
	win[1].Double(Q)
	for i := 3; i <= 15; i += 2 {
		win[i - 1].Add(&win[i - 2], Q)
		win[i].Double(&win[((i + 1) >> 1) - 1])
	}

	// Initialize P with the top digits (which are all nonnegative).
	j := int(sd0[25])
	if j == 0 {
		*P = jq255sNeutral
	} else {
		*P = win[j - 1]
	}
	P.addFromWindowAffineVartime(&jq255sWin_G_eut, sd1[25])
	P.addFromWindowAffineVartime(&jq255sWin_G130_eut, sd1[51])

	// Process all other digits in top to bottom order.
	for i := 24; i >= 0; i -- {
		P.DoubleX(P, 5)
		P.addFromWindowVartime(&win, sd0[i])
		P.addFromWindowAffineVartime(&jq255sWin_G_eut, sd1[i])
		P.addFromWindowAffineVartime(&jq255sWin_G130_eut, sd1[26 + i])
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
func (P *Jq255sPoint) MapBytes(bb []byte) *Jq255sPoint {
	// Decode source bytes as a field element. This applies modular
	// reduction.
	var f gf.GF255s
	f.DecodeReduce(bb)

	// Formulas: we map into a point on the dual curve
	// y^2 = x*(x^2 + aa*x + bb) with aa = -2*a and bb = a^2 - 4*b.
	// We use a fixed non-square constant d.
	//
	// Then, the two candidates for x are:
	//   x1 = -aa / (1 + d*f^2)
	//   x2 = -x1 - aa
	//      = -aa*d*f^2 / (1 + d*f^2)
	// The candidates for y^2 are then yy1num/(yden^2) and
	// yy2num/(yden^2), with:
	//   yy1num = -aa*bb*d^3*f^6 + (aa^3-3*aa*bb)*d^2*f^4
	//            + (aa^3-3*aa*bb)*d*f^2 - aa*bb
	//   yy2num = -aa*bb*d^4*f^8 + (aa^3-3*aa*bb)*d^3*f^6
	//            + (aa^3-3*aa*bb)*d^2*f^4 - aa*bb*d*f^2
	//          = yy1num*d*f^2
	//   yden = (1 + d*f^2)^2
	// If yy1num is a square, then we use y = sqrt(yy1num) / yden.
	// Otherwise, we use y = -sqrt(yy2num) / yden.
	//
	// Once (x,y) is obtained, we apply the theta_{1/2} isogeny to
	// get a point in the right group.
	//
	// If 1 + d*f^2 == 0, then we map to the neutral N.
	//
	// For jq255s:
	//   a = -1   b = 1/2   aa = 2   bb = -1   d = -1
	//   x1num  = -2
	//   x2num  = 2*f^2
	//   xden   = 1 - f^2
	//   yy1num = -2*f^6 + 14*f^4 - 14*f^2 + 2
	//   yy2num =  2*f^8 - 14*f^6 + 14*f^4 - 2*f^2
	//   yden   = (1 - f^2)^2
	// Since we want u = x/y, we have:
	//   u = U/T = (xnum*yden) / (xden*ynum)
	//           = (xnum*xden) / ynum
	// since yden = xden^2.

	var f2, f4, f6 gf.GF255s
	f2.Sqr(&f)
	f4.Sqr(&f2)
	f6.Mul(&f2, &f4)

	// yy1num and yy2num
	var yy1num, yy2num, tt1, tt2 gf.GF255s
	tt1.Sub(&f4, &f2)
	tt2.Lsh(&tt1, 3).Sub(&tt2, &tt1)
	tt1.Sub(&tt2, &f6).Add(&tt1, &gf.GF255s_ONE)
	yy1num.Lsh(&tt1, 1)
	yy2num.Mul(&yy1num, &f2).Neg(&yy2num)

	// Assemble the point coordinates; this is still on the dual
	// curve. If f != +/- 1, then it cannot happen that yy1num == 0.
	// However, it is possible that yy2num == 0, but only if f == 0;
	// in that case, yy1num == 2, which is not a quadratic residue,
	// and yy2num is selected.
	ls := yy1num.Legendre()
	qr1 := 1 - (ls >> 63)
	var xnum, xden, unum, uden gf.GF255s
	xnum.Neg(&gf.GF255s_ONE).Select(&xnum, &f2, qr1).Lsh(&xnum, 1)
	xden.Sub(&gf.GF255s_ONE, &f2)
	unum.Mul(&xnum, &xden)
	uden.Select(&yy1num, &yy2num, qr1).Sqrt(&uden)
	uden.CondNeg(&uden, 1 - qr1)

	// If f == +/- 1, then (xnum:xden:unum:uden) == (+/-2:0:0:0)
	// If f == 0, then (xnum:xden:unum:uden) == (0:1:0:0)
	// Otherwise, xden != 0 and uden != 0, and the value is correct.
	to_fix := uden.IsZero()

	// Apply the isogeny theta_{1/2} to get a point in the proper
	// group on the right curve.
	//   x' = 4*b*u^2
	//   u' = 2*x/(u*(x^2 - bb))
	var xn, xd, un, ud, en, ed gf.GF255s
	xn.Sqr(&unum).Lsh(&xn, 1)
	xd.Sqr(&uden)
	un.Lsh(&uden, 1)
	xden.Sqr(&xden)
	ud.Sqr(&xnum).Add(&ud, &xden)

	// If we are in a special case (f = 0, 1 or -1), then un = 0 and
	// ud != 0 (which is correct), but xn and xd are incorrect and
	// we want to fix them.
	xn.Select(&gf.GF255s_ZERO, &xn, to_fix)
	xd.Select(&gf.GF255s_ONE, &xd, to_fix)

	// Convert to EZUT coordinates.
	var k1, k2 gf.GF255s
	k1.Lsh(&xn, 1).Sub(&k1, &xd).Mul(&k1, &xn)
	k2.Sub(&xn, &xd).Mul(&k2, &xd)
	en.Add(&k1, &k2)
	ed.Sub(&k1, &k2)
	var ud2, uned gf.GF255s
	ud2.Sqr(&ud)
	uned.Mul(&un, &ed)
	P.e.Mul(&en, &ud2)
	P.z.Mul(&ed, &ud2)
	P.u.Mul(&ud, &uned)
	P.t.Mul(&un, &uned)

	return P
}

// =====================================================================
// Precomputed windows in extended affine (e, u, t) coordinates:
//   i*G         for i = 1..16
//   i*2^65*G    for i = 1..16
//   i*2^130*G   for i = 1..16
//   i*2^195*G   for i = 1..16

// Points i*G for i = 1 to 16, extended affine coordinates
var jq255sWin_G_eut = [16]jq255sPointAffine {
	// G * 1
	{
		e: gf.GF255s { 0x104220CDA2789410, 0x6D7386B2348CC437,
		               0x55E452A64612D10E, 0x0F520B1BA747ADAC },
		u: gf.GF255s { 0x0000000000000003, 0x0000000000000000,
		               0x0000000000000000, 0x0000000000000000 },
		t: gf.GF255s { 0x0000000000000009, 0x0000000000000000,
		               0x0000000000000000, 0x0000000000000000 },
	},
	// G * 2
	{
		e: gf.GF255s { 0x155215A00E9EA08D, 0xC474598BF46D0A34,
		               0xA479391BDE7F0296, 0x7DCAB2C9EFDB7348 },
		u: gf.GF255s { 0xB3E22F8D0D1657FC, 0xE5427944219E490E,
		               0x256C2BE75887FD30, 0x6F44EC749C5869ED },
		t: gf.GF255s { 0x84CC11AA69B07916, 0xEC33C75916FEEF18,
		               0x2501ACD978762D61, 0x42B401D3D5F7C6BD },
	},
	// G * 3
	{
		e: gf.GF255s { 0x1C2AE84D369F681F, 0xE7931687A657BB05,
		               0xFAC101D18780329D, 0x407A13E6346F39AA },
		u: gf.GF255s { 0x7204233F36F06441, 0x36F07204233F36F0,
		               0x233F36F07204233F, 0x7204233F36F07204 },
		t: gf.GF255s { 0xB582FD5064ED4BDD, 0x9F6417AAAE3D0228,
		               0x10163C9F79D25F29, 0x3E21C102D9DADBFA },
	},
	// G * 4
	{
		e: gf.GF255s { 0xFCB9E72D0A0C8EFE, 0x879B715FCBC752B5,
		               0xADBF0638585A38A4, 0x1153E6DF83525815 },
		u: gf.GF255s { 0x9204A59E69223E39, 0x4E645F874B12D8E7,
		               0xF2A1145C9F5D3475, 0x54E64904662F1657 },
		t: gf.GF255s { 0x330B2867F87650FD, 0x013E9EBA0D1CB6FA,
		               0xAB4BE93A7FC59C00, 0x11F4E2818B86CA68 },
	},
	// G * 5
	{
		e: gf.GF255s { 0xDE7D4909B7768768, 0x89E1BB43C57ADC57,
		               0x850ECBC936A0B118, 0x67D128B9D764C85E },
		u: gf.GF255s { 0xDF0337C00667B64D, 0x58856B292FBA673A,
		               0xC17C3E9333A6D7CE, 0x52932B9A9F0CC65D },
		t: gf.GF255s { 0x1C824AAF423C38C3, 0x635DC68819778FD9,
		               0x3B8EBB2DDCBD289A, 0x591738A93D12426F },
	},
	// G * 6
	{
		e: gf.GF255s { 0xDB908132D0DC5ED9, 0xBAF1C9AE563AB60C,
		               0x3F0A5A3647C0974E, 0x1AE40CE8E4758D2C },
		u: gf.GF255s { 0x2378FCE7659F8304, 0xF71199A7E92DA598,
		               0x05D59CECEE1D7E76, 0x7B55CE1D18E39726 },
		t: gf.GF255s { 0x5C964B33DBE905C9, 0xFB8008F322A48DBC,
		               0x002552F2D3C7EDD4, 0x42B373007F5A5422 },
	},
	// G * 7
	{
		e: gf.GF255s { 0xBCF616435A305CD1, 0x8FF2EAAA7A24E5FE,
		               0x1F836A920FFADF5D, 0x778691038A25A05C },
		u: gf.GF255s { 0xBB70A3099712F248, 0x62AE8CABB5C7CED6,
		               0x150E35D2C3D060D0, 0x6EB76B2647D95DA4 },
		t: gf.GF255s { 0x5429F9F2E42EC752, 0x64F44E4B72EDC8DB,
		               0x6A2EACE6BC1B7A38, 0x59A0770044734AB4 },
	},
	// G * 8
	{
		e: gf.GF255s { 0x1EA227FB56AF9D20, 0x342B4BA3A435E2AC,
		               0x188E0332E4A5FB16, 0x4CB7E380E64F8650 },
		u: gf.GF255s { 0x2DE04D6F93F6EEA0, 0x13A2A2DE4BD9AB93,
		               0x4BA8485FACF9CC03, 0x26DCD74FEC331ED1 },
		t: gf.GF255s { 0x008F881CFF3A091D, 0xD79E9196AC529912,
		               0x27641FC9B91D9EE2, 0x08F2F46F791D0D36 },
	},
	// G * 9
	{
		e: gf.GF255s { 0x7064AA193867E249, 0x9E0262DD1BBD9331,
		               0x28DAB54C47465363, 0x186E71D5F58996EC },
		u: gf.GF255s { 0xF293B01428B2AEF7, 0x6544BF679F64AADF,
		               0xC25DB7DB3DC0038B, 0x641BC3EAE8B16348 },
		t: gf.GF255s { 0xC9769297553D878E, 0x5AFEBCBFF210B104,
		               0xC30E7E91EBC5E399, 0x5344F6CBF78889DE },
	},
	// G * 10
	{
		e: gf.GF255s { 0xEC41FB9869A38BAF, 0x55781708B257A28D,
		               0xF1BEAC3346435A2C, 0x24D9FE2BF2E5C1DF },
		u: gf.GF255s { 0xAC2450CD35E0E33E, 0xA5FBAFCA4BE098FF,
		               0x15C5CFFBF6CC3634, 0x11EB002E4DD0CC8A },
		t: gf.GF255s { 0xB57551F40B3E9878, 0xA1D9E1AC087CD2FE,
		               0xFE443CA5E072548A, 0x3361CE26431C24BB },
	},
	// G * 11
	{
		e: gf.GF255s { 0x640788E8FDB68B79, 0xA02A4D6FD4645856,
		               0xB18CAB0F6A54EB86, 0x01511890814FCFF1 },
		u: gf.GF255s { 0x4845EE2F9FDD13DE, 0x343EAC1284D73CA8,
		               0xD7761D7533C3A9ED, 0x44204FCFD6AF1B44 },
		t: gf.GF255s { 0x88C209C7EA04605B, 0xDEC8E3E4B025867C,
		               0xF25BD01AAB8AE0C4, 0x10E31B1DBFBE67E2 },
	},
	// G * 12
	{
		e: gf.GF255s { 0xC31D48EEAC5599AB, 0x3E87D66869F1BBA3,
		               0x8AB438671165365F, 0x42B172782E25AB56 },
		u: gf.GF255s { 0xEDFD8B4BAD8160D1, 0x9BBB9D365457C4A7,
		               0xC2A3F9623BA7F9DB, 0x1EA722FF296F8D92 },
		t: gf.GF255s { 0xDBD377FD4727F1F8, 0x060E2D9DB05F1736,
		               0xA86F76D59422C162, 0x466ED47592E58063 },
	},
	// G * 13
	{
		e: gf.GF255s { 0xFFAD311FC42DF019, 0xCC5F0E316B1677A5,
		               0xDB27E2AA724C430C, 0x72491B477910438F },
		u: gf.GF255s { 0x4DEC233876049D71, 0x12E5410DA8EA7661,
		               0xA8C80C9E1E78A79B, 0x56187141071FC1D5 },
		t: gf.GF255s { 0x5F87D7B161AB363D, 0xFADED89F6751C8BA,
		               0x37A78C21BA3C4BCD, 0x3E2BE3EF9C950F96 },
	},
	// G * 14
	{
		e: gf.GF255s { 0xB1620E564FE98F30, 0xB119CC9192E72DD0,
		               0xA0040674E6CCAC7F, 0x5E430A47DD41A2D1 },
		u: gf.GF255s { 0x0693069387390957, 0xC3E0CD34567E3C9A,
		               0x2EDB82538E0CEEAE, 0x60EF0D216D3123EB },
		t: gf.GF255s { 0xCEB292E0EC603124, 0x18B961A309C8339D,
		               0x917BC5415D3B524A, 0x681E3041257F38AE },
	},
	// G * 15
	{
		e: gf.GF255s { 0x871FFDDBCBB41565, 0x98B526861DE23B6B,
		               0x51EBFEAAC016B5CD, 0x3BE6B06EB08C48A2 },
		u: gf.GF255s { 0x3010093957EA5D5A, 0x99E3DE25AE95E20B,
		               0x5594A4E0449B9EFD, 0x100615BCCFB9A5A5 },
		t: gf.GF255s { 0x2CB15E4AD534F12E, 0x07A2C5FB49523FF8,
		               0x79BCF0468E59F7E5, 0x172AAB628A8B99AA },
	},
	// G * 16
	{
		e: gf.GF255s { 0x18BF1FBA94943C22, 0xEE888DF3C9AF901F,
		               0xE6355FA5083B2FE3, 0x2DE21A5F5B5F3676 },
		u: gf.GF255s { 0xFE9991C8BEE08226, 0x70631241B6132ED0,
		               0xF17032C1930CC3DF, 0x1880F6E65E61EDCA },
		t: gf.GF255s { 0x7D5A50285053E953, 0xB0627DFEE2C9B21E,
		               0xB6548C27E54A4C3C, 0x5C75B961A2BDDB5E },
	},
}

// Points i*(2^65)*G for i = 1 to 16, affine extended format
var jq255sWin_G65_eut = [16]jq255sPointAffine {
	// (2^65)*G * 1
	{
		e: gf.GF255s { 0xAD671792FC850F1F, 0xB7D10E9710F4FB5A,
		               0xE59F700CF1498B64, 0x145C3F79A555FEFD },
		u: gf.GF255s { 0xC504DF70301870EC, 0xC3F5757559BEB30B,
		               0x59EF9CAA0E041627, 0x3AA6C1241FDF29EA },
		t: gf.GF255s { 0x6D1B6960A73A5D82, 0x7E905634B7B83754,
		               0x5728E3A866E7AD73, 0x6F9938F04DD20707 },
	},
	// (2^65)*G * 2
	{
		e: gf.GF255s { 0xAEA44289F769D548, 0xF3B3A271D9C85032,
		               0x450F6EC619A10A0D, 0x2A078A32D8A60452 },
		u: gf.GF255s { 0xE395C342A49AB090, 0xBEEAEA12D994BBBA,
		               0xE555C4B2316756DE, 0x0C4B1FDF2A2835EA },
		t: gf.GF255s { 0xACC1B6C0E71A2679, 0x3B828FBF6A20F966,
		               0x2B40BE24F305035E, 0x5138EFFC5F9EE71D },
	},
	// (2^65)*G * 3
	{
		e: gf.GF255s { 0xFCB8737061C3B389, 0x2EB5074C5FB663E6,
		               0xC58DE35B727B63A1, 0x6F9AECFE00B2FBF9 },
		u: gf.GF255s { 0xDD2736E064D47F35, 0x8235DE8C4CCCB7AD,
		               0x84993635937FCD8B, 0x05DADA9DFEACCC60 },
		t: gf.GF255s { 0x041CF837CCDEEA2C, 0x128FA18D8F430D4E,
		               0x146AA2BAF6EC5EF7, 0x367779B57A3314D7 },
	},
	// (2^65)*G * 4
	{
		e: gf.GF255s { 0x02E353146AB767DD, 0x07B12ACEF6DB6D43,
		               0xB844A3A38291AF02, 0x168FE09EBFF84166 },
		u: gf.GF255s { 0x20CCC8942E0DD4AF, 0xB05A144640E6146F,
		               0x39FB3501747B5584, 0x637DD68BAA0DCC5D },
		t: gf.GF255s { 0xF858BE49C8D19D78, 0x1CD2F537C049D501,
		               0x5B82CA8E72C504DC, 0x081DAA0CD8F74B16 },
	},
	// (2^65)*G * 5
	{
		e: gf.GF255s { 0xEAAC278EE492D81E, 0x9B679600DF4B5715,
		               0x0930AFE6F7D251E0, 0x12B18E684B0A57E8 },
		u: gf.GF255s { 0x8E3257DDE0B61BA4, 0x034912BE87333A80,
		               0x355389746FE3860D, 0x0430DD4DB72FDAE3 },
		t: gf.GF255s { 0x8F69AB4C10703FE2, 0x9DB06D02B9CA5DEA,
		               0x80E1FCACE9A31DF3, 0x6CE5827F833EE541 },
	},
	// (2^65)*G * 6
	{
		e: gf.GF255s { 0x296B5AC7D0BE8BD1, 0xF9ABC6FE17A925CB,
		               0x9A382C9BFE2755E0, 0x18BD1D68A45DEDDC },
		u: gf.GF255s { 0x41C6037317995DE4, 0xEB9F7DAEFD0E5177,
		               0x3C89F25AEBEFEA23, 0x6BC51220E2A6FEA1 },
		t: gf.GF255s { 0xA8BE5154DC0850D8, 0x46B4997F9117FF61,
		               0xF3632B4CFD05C7A5, 0x5198473D3A06C5F9 },
	},
	// (2^65)*G * 7
	{
		e: gf.GF255s { 0x0C025505E54EC194, 0x16DFAAD15E54A341,
		               0x4E0E09389BAC268D, 0x2C841A6C3AE0F758 },
		u: gf.GF255s { 0x6E0C6EC00B45CB36, 0xDE22614DA09D7B7D,
		               0x68532DB386C5311B, 0x6C76A36690BF3721 },
		t: gf.GF255s { 0xEF4DC9C561021E16, 0xBF20C5FEEE045148,
		               0xCC0911EBF7E93EA0, 0x2C44B584B042AE8F },
	},
	// (2^65)*G * 8
	{
		e: gf.GF255s { 0x0D409883C7DFAAFF, 0xA9BB3C1344A2B7AC,
		               0xADEE45CF23AA330B, 0x32F8215138DF8DF9 },
		u: gf.GF255s { 0xCF5B38B898027EEE, 0x2E0E5F286C358DFB,
		               0x77F9F481B54C83A0, 0x320289A4D4AD38D8 },
		t: gf.GF255s { 0x850595A42D050368, 0x69338F651A87EED4,
		               0xB56B41AAC63C7E31, 0x0A7DECD4F28D875D },
	},
	// (2^65)*G * 9
	{
		e: gf.GF255s { 0x3EA9CA530CD6C445, 0xF4121ED9B11E2F78,
		               0x32D2ACB4DF4CBE2B, 0x6BFAD4535C93F764 },
		u: gf.GF255s { 0xE9F73ACE30D70239, 0x2D6EB28373FFD2AB,
		               0xE362CB1B14264D84, 0x6A58C4C447F58EDA },
		t: gf.GF255s { 0xD421A3CDEB81C4F0, 0x8C5750AC973B8F0F,
		               0x43D7F20CB3184DFA, 0x05EA57EAA00332B5 },
	},
	// (2^65)*G * 10
	{
		e: gf.GF255s { 0x6E80B93F8C55EFC5, 0xCFA05077AC82CA7F,
		               0xFC3A41A8997B623A, 0x594955BFF388397E },
		u: gf.GF255s { 0xBD4ABD3A1D23F58F, 0xFDF2BE7E13D824B7,
		               0xFDF4EBAD6DA9D7C5, 0x6F6BF0FD0A760CE6 },
		t: gf.GF255s { 0xB81C256655EADF3C, 0xBDD6F1B6D18BC1A6,
		               0x3160D36F471E8803, 0x5B172E38809C36C6 },
	},
	// (2^65)*G * 11
	{
		e: gf.GF255s { 0x80DDD82DF53D8310, 0x9FE3E489FB0E0E43,
		               0xE1EC96F492FEAA7B, 0x0D9EFEC28CA0A54C },
		u: gf.GF255s { 0xAC82BC835BEF5D82, 0x1ED56E4F26768EEA,
		               0x7CB1DF78AE0F6520, 0x3F27081B1A793D69 },
		t: gf.GF255s { 0x057471C1F93EB90E, 0xDA034A47FB299F2D,
		               0xC51142A5B0691BD2, 0x5FD614C1FCF6A3EE },
	},
	// (2^65)*G * 12
	{
		e: gf.GF255s { 0xDBEFD63B52CB2A5C, 0xCA182B05C7FCE11E,
		               0x8173D517B7E0B2B0, 0x1A664CA56EBCC5BD },
		u: gf.GF255s { 0x3C2CD45E00488C1B, 0x744209C1D5AA1E98,
		               0xB0A6FCE83F628FFF, 0x770D858C1138F9B1 },
		t: gf.GF255s { 0x6A96C56588D37A03, 0x733C7CA714153EC3,
		               0xA04AF9D178D61F3C, 0x10420882FFE4F21E },
	},
	// (2^65)*G * 13
	{
		e: gf.GF255s { 0x855811311CCCC730, 0x7B6D5AAC61E9960D,
		               0xC3729090E26F9B9C, 0x035DE04A61BABE1E },
		u: gf.GF255s { 0x2B0FECC1159B7C7E, 0x5176B2243AD75522,
		               0xF578275E6C3BD8ED, 0x3EA45CDA4A4335B5 },
		t: gf.GF255s { 0x399CEF5FFB8296AD, 0x7BB775A10203F23D,
		               0xEA10679B009BFFEE, 0x70977639FA69E967 },
	},
	// (2^65)*G * 14
	{
		e: gf.GF255s { 0x969B92FA2E20E4A3, 0xC10608CEDBFB267E,
		               0x03B072171B964F50, 0x513E12D64D4EAAF8 },
		u: gf.GF255s { 0xDCABCD1C2B792081, 0xC0B470A297317B12,
		               0xB2C96A085D26A03D, 0x52A2BF2E7D7202E9 },
		t: gf.GF255s { 0x92EA95A74111A128, 0x8164FFF4195F29D9,
		               0x8A3F159339E5051E, 0x43FC22627E870D74 },
	},
	// (2^65)*G * 15
	{
		e: gf.GF255s { 0x5B4B80D1C3DF479B, 0x26E5E74F1817958B,
		               0x899B72B33AED803B, 0x07D749F9FAB7DFA1 },
		u: gf.GF255s { 0xCB4E8E8682ED1F56, 0x95E7B3DAE6F438C6,
		               0x436574C311E4EAE5, 0x5B29A9D46B6557D1 },
		t: gf.GF255s { 0xA0FF9F128E2A8D1F, 0xD3C936A9527B4E88,
		               0xB279CEBEF2B0454A, 0x04B77CA16130D747 },
	},
	// (2^65)*G * 16
	{
		e: gf.GF255s { 0xD82E7D79DB9A5B83, 0xE55F1507DB1386D1,
		               0x4E2ED8BC0B19751E, 0x191D614E6AD8C891 },
		u: gf.GF255s { 0xA0A662EF07947CAE, 0xCFD4B9B4407ADF90,
		               0x9698A42E4CDA26D4, 0x00D98818FA1C5F43 },
		t: gf.GF255s { 0xF048C3E420D474E3, 0xA09CE0B0CDA1B6AB,
		               0x31ED70E89727A4CA, 0x214C5B3CB2D89DF0 },
	},
}

// Points i*(2^130)*G for i = 1 to 16, affine extended format
var jq255sWin_G130_eut = [16]jq255sPointAffine {
	// (2^130)*G * 1
	{
		e: gf.GF255s { 0xC547E3D2287B7A8C, 0x907B64252E9A0A54,
		               0x7A7BA22355F5D398, 0x5FEB95EE492D6D08 },
		u: gf.GF255s { 0xC17AB82D247C18A0, 0x95542A3E6973F13C,
		               0xB14CDFC79E957BD2, 0x661229C7BADE4F32 },
		t: gf.GF255s { 0x5541C59928E441B7, 0x5B912D60442F368E,
		               0x2EC475815BC7D073, 0x02802A6A03E1A1E7 },
	},
	// (2^130)*G * 2
	{
		e: gf.GF255s { 0x15A644D90F28E7B4, 0x2AEBD6074550787C,
		               0x319418ED83354F29, 0x2B3CECA02688A0B0 },
		u: gf.GF255s { 0x72C27DF187129578, 0xF443CA1C94A3CAEB,
		               0xEA24368F35C3A22D, 0x17A619CA7283DADF },
		t: gf.GF255s { 0x0EC06DB43459B171, 0x1A6558F98FD7DD2D,
		               0x4B3F9D68A82A6DEE, 0x649637A5340F9542 },
	},
	// (2^130)*G * 3
	{
		e: gf.GF255s { 0xA065A1038F86DA38, 0xBF1723E466338796,
		               0xD922BF8C18633561, 0x168AAAD797CE7EDA },
		u: gf.GF255s { 0x31D729A12CB400ED, 0x7A195520CEEAA6A9,
		               0x86F51D53192BE4B6, 0x22A23BF61CDCF306 },
		t: gf.GF255s { 0xCACA700AE212FCC2, 0x8E87630B3381880A,
		               0xC73BF67BBC73A980, 0x66C538B71AA1CCE1 },
	},
	// (2^130)*G * 4
	{
		e: gf.GF255s { 0xA7C985807B396AC1, 0xF5B27B3B050B1CC0,
		               0x199AF9DFF3BB8550, 0x0DD0B15BCCECE2DE },
		u: gf.GF255s { 0x767FEA8208A28FA3, 0xA2DC84A037DCA7EE,
		               0xEBD2AC8233E1DE6F, 0x188F708521AD21CA },
		t: gf.GF255s { 0x3AC52264BA423FFD, 0xEF5048FF3B333AA4,
		               0x2FBD912DEC2F8A1A, 0x2BEF1130670533CF },
	},
	// (2^130)*G * 5
	{
		e: gf.GF255s { 0xCD8A44DA9FC7024F, 0xF737A550BC29AA93,
		               0x0D8952811E9106B9, 0x2392F5F42DA132F4 },
		u: gf.GF255s { 0xB820181764875239, 0x5B4DD1BB342A7459,
		               0xC66B48D80BE03B56, 0x7F08F18F636EA208 },
		t: gf.GF255s { 0x2CC6020672019CF2, 0x14E86A7605B43C11,
		               0xB9CE10F55DE4C795, 0x67B29276BA073D7D },
	},
	// (2^130)*G * 6
	{
		e: gf.GF255s { 0x43D22D5EF864F06F, 0x69548FFF3BE07C7C,
		               0xEBB4F5A33EDD77B8, 0x5C5A93242CA0CD5B },
		u: gf.GF255s { 0x9FFA28B92E2ECC30, 0xFB480B5678DA5C00,
		               0xE1C94F9AB14804D7, 0x7D453CC948C9FA1D },
		t: gf.GF255s { 0x99B10D2684BD090A, 0xAA8361C9702948B9,
		               0x8FE833DB15F7685C, 0x5A40858A46BA5349 },
	},
	// (2^130)*G * 7
	{
		e: gf.GF255s { 0xFE714DB9C806064F, 0xEC532CBC96052BB3,
		               0xFEB5667662B03E4E, 0x7D201AFC52734277 },
		u: gf.GF255s { 0x8784811955F5DEF5, 0x48D2662133271AC2,
		               0xAF3F2A3FA6C81237, 0x1A086A6FE01B2E7D },
		t: gf.GF255s { 0xE60717C5AA0289C2, 0x06226B3A465F8A77,
		               0xAAE6A9C9C053632B, 0x7444100416332425 },
	},
	// (2^130)*G * 8
	{
		e: gf.GF255s { 0x9C901178C3D8E685, 0x3FD38D2D5E9E8C93,
		               0xA0478BE33AECC482, 0x507DF5A8263956AC },
		u: gf.GF255s { 0x1DC9A2A859DBD92B, 0xBF15C84C6245BBB3,
		               0x2A8D65D0B2CAB143, 0x3DE43A8017A5525E },
		t: gf.GF255s { 0x2040089013E8DD16, 0x98B672949E66D5ED,
		               0x2468A4B876BC1D61, 0x1D2C24FA4A76A7C6 },
	},
	// (2^130)*G * 9
	{
		e: gf.GF255s { 0x4E0B3042A15BFEC5, 0x25107B05B78EFCC3,
		               0x895F51290AB84B6A, 0x2875F4526CD0570A },
		u: gf.GF255s { 0xDB52105D891C8CF1, 0xA1E811D5F01D372F,
		               0x5CF867DE0ADED951, 0x17052800E0B4DA8B },
		t: gf.GF255s { 0x465AF151D990BDFE, 0x61C1C858ED9253B0,
		               0x565E0443A3FBDE66, 0x12DB650B53E17835 },
	},
	// (2^130)*G * 10
	{
		e: gf.GF255s { 0x0149667229D2E833, 0x01916903ED0CA017,
		               0xB4E997E9293857A0, 0x52B3D9F485A3DB38 },
		u: gf.GF255s { 0x126C847746496B31, 0x23D37316626487D7,
		               0x46505E02982B39EB, 0x234345752143582E },
		t: gf.GF255s { 0x6BBD04A3A0420BAA, 0x9684CAEB425C1B2E,
		               0x2D86F80589C6CAF4, 0x625467C6A6568CBE },
	},
	// (2^130)*G * 11
	{
		e: gf.GF255s { 0x8444EE09228C2235, 0xC67ADEEDBBF1E2D0,
		               0x78F69BE716B2CFD7, 0x7B0AF5F2B95FF650 },
		u: gf.GF255s { 0xD51EDD52BECE2651, 0x792F3E2ED6FA4957,
		               0x8349ED268E6F750E, 0x6BB94CDFEF6F6D88 },
		t: gf.GF255s { 0xE4774F545088C5A8, 0x43CAA7ECFFA8E04D,
		               0x2CBAD6B90644B58D, 0x66EC80D6A1180C6F },
	},
	// (2^130)*G * 12
	{
		e: gf.GF255s { 0xAF4B65BA8591C1AD, 0x3494662617CEAF23,
		               0x67D0F7B0D86050C8, 0x7DD5EFEF2F7E3A00 },
		u: gf.GF255s { 0x89340045DE2129A6, 0x10CF7923E46000DB,
		               0xA17D1BEB2C69ABE2, 0x1F3592DE203C2DA4 },
		t: gf.GF255s { 0x946A6220C9A50DDD, 0xB19B06C590529CD9,
		               0x932CF781C2D4D2BA, 0x643EC0482633F97D },
	},
	// (2^130)*G * 13
	{
		e: gf.GF255s { 0x681F8DB202E72852, 0x1321374B6FBEA675,
		               0xD874B17723DB7E8C, 0x19430AF80693A8A2 },
		u: gf.GF255s { 0x0F4767AD57355C90, 0xCAD900D819124EA8,
		               0xB4B045E5E702318A, 0x7AB6CB353A7E1058 },
		t: gf.GF255s { 0xF5572A210505358F, 0x0C528F0B20133EC8,
		               0xC425B464B3606AA3, 0x33E0CFF16CAB14B8 },
	},
	// (2^130)*G * 14
	{
		e: gf.GF255s { 0x1739D8DE862E3AF0, 0xEAEA8FAEB3D5ED7B,
		               0x394668C34B06CE7A, 0x26B43F065FC95D4F },
		u: gf.GF255s { 0xE62B26D301A95CF8, 0x84BC9EEC48398A95,
		               0x08F76F6E267875A4, 0x32D3B1A49A0A50A1 },
		t: gf.GF255s { 0x1E1D0050EAB69B8C, 0x5FAE0C5B31875B73,
		               0x5FD8D77B29295DB8, 0x720DA33C8F5136E0 },
	},
	// (2^130)*G * 15
	{
		e: gf.GF255s { 0x53614A3BC609C455, 0x1FDF44E6D229E43B,
		               0xC69D5CC6EC35FF3A, 0x40D8B2EDD4386D2C },
		u: gf.GF255s { 0x9A0B63FF58819185, 0xF789BF221086D125,
		               0xAD341ED8B1E1776A, 0x082617E146EBF733 },
		t: gf.GF255s { 0x94BAFDC6B0AF72AF, 0xCE7FD95867032F80,
		               0x60A041C742E3E324, 0x52BAFECDE0C48A92 },
	},
	// (2^130)*G * 16
	{
		e: gf.GF255s { 0xB59722476686F776, 0xFCAABBE442C88D93,
		               0x0A7A9E0DFA237DD5, 0x2DBB938F43FCED53 },
		u: gf.GF255s { 0xF50D8973C1409964, 0x7E7D988FECBBDB1D,
		               0x64BAB1680DE07B83, 0x588C179EE277A32E },
		t: gf.GF255s { 0x62A6B63BF1FFD3B3, 0xAA7B23F4A4AEC308,
		               0x52196FBDBC75755C, 0x60C42361D4635289 },
	},
}

// Points i*(2^195)*G for i = 1 to 16, affine extended format
var jq255sWin_G195_eut = [16]jq255sPointAffine {
	// (2^195)*G * 1
	{
		e: gf.GF255s { 0xEB4B028B52EF088D, 0x7BFFA172511D0BF8,
		               0x20A728914AA70C44, 0x1982687F23A0A46D },
		u: gf.GF255s { 0xD5386C024F6950E1, 0xB46F6A19F2C711D7,
		               0xD13612F77F90F84F, 0x02B1FE40C35B1108 },
		t: gf.GF255s { 0x24BC6E863C8937C9, 0xA15165FBD4AE4DC9,
		               0xF54207D14BAA772F, 0x4ECAC7AF734C65E6 },
	},
	// (2^195)*G * 2
	{
		e: gf.GF255s { 0x6736CB36F5DF3C6B, 0xAAC155C49E622FFF,
		               0x08BF6AEC6C9AF425, 0x32453BFF3E754687 },
		u: gf.GF255s { 0x4D20522801D9BE69, 0x491F6D9530005684,
		               0x0394879561348B5B, 0x207807A3C7DA4C9E },
		t: gf.GF255s { 0x39C460145081D09A, 0xE3BB41D1EEDD6F53,
		               0x8BAB5AC42AECE484, 0x7BE8AD2F073E261D },
	},
	// (2^195)*G * 3
	{
		e: gf.GF255s { 0xE2163ED40240D5AB, 0x0D36DDAFD2EDE213,
		               0xB5DA9A3A84DA7A94, 0x616BB17981AD5800 },
		u: gf.GF255s { 0x62329BF642399CCB, 0xD15830A754BDD8E5,
		               0x052C83B745F76715, 0x2C2FEC4277A075E5 },
		t: gf.GF255s { 0xCA181F67050E4AFB, 0x8FA110DC3CABC581,
		               0xEC7B337F46B0A4C3, 0x614659258F72400F },
	},
	// (2^195)*G * 4
	{
		e: gf.GF255s { 0x11A5ECC6C01F036B, 0x00905EC73DCFAB11,
		               0xBD182106BE833EB9, 0x21345004520D6147 },
		u: gf.GF255s { 0x33F1C936823F9793, 0x45EB4CD78646E9E2,
		               0xB28C3C8D49110514, 0x20B293F739F3AA65 },
		t: gf.GF255s { 0xECE0E7C9F465E491, 0x6E3D99DD3E295315,
		               0x11CAE9FC0D23AB0D, 0x395BB5CD3D65F08E },
	},
	// (2^195)*G * 5
	{
		e: gf.GF255s { 0xC3465879B04EFBCB, 0xDDD7400D7DF1D9F8,
		               0x3FC45AA3BB68E927, 0x0D654BC97CC560D9 },
		u: gf.GF255s { 0x3F09E8FE4D4F0DAA, 0xE47EE2E2B553CA96,
		               0x9FFA3C3AAC578533, 0x356B0675A08113C1 },
		t: gf.GF255s { 0x2F2E5ED16CF3E78E, 0x3739979BDD5F3D89,
		               0x6DA80DAB9FB489BE, 0x431D2978BDA14A65 },
	},
	// (2^195)*G * 6
	{
		e: gf.GF255s { 0x93683AAAE71501A7, 0x7A150F1801C67235,
		               0x5F345C36300E35B7, 0x7A900684605D8B43 },
		u: gf.GF255s { 0xA143C257C8E7B387, 0x62D889136CC108D4,
		               0xE32F52E0BB72837C, 0x076CCDEE12F67001 },
		t: gf.GF255s { 0x02F43BF7262D9E46, 0xC600EF045B679BDB,
		               0x75B858E15A115558, 0x580DEEF025934593 },
	},
	// (2^195)*G * 7
	{
		e: gf.GF255s { 0xCAF3D56D09001CAE, 0x56CCAFC57F61B1E3,
		               0x6BE0EB313BA1A346, 0x4AB3F0BEB592478D },
		u: gf.GF255s { 0x971F761E8AEDB101, 0xD10EB0D73B802862,
		               0x100B543BDE5AD3F6, 0x2D2D355F27445B8E },
		t: gf.GF255s { 0xBB6B514090669B3D, 0x6AF732299368BB18,
		               0xD33B393FC4662114, 0x44B393C42E5F9548 },
	},
	// (2^195)*G * 8
	{
		e: gf.GF255s { 0x073F794DC96EBE06, 0xF9AAA7896642E5AB,
		               0x47F0972818DC4F59, 0x30E51AB070B031C5 },
		u: gf.GF255s { 0x6AC9953979CB6EA2, 0xB754075F43904E3A,
		               0x0F5C5D706C51FC54, 0x72F0811D6A945237 },
		t: gf.GF255s { 0xBB17DDD37656E362, 0x9F321433405FA3BF,
		               0x326CCE6372A57EE6, 0x2CEB4E800F2A9875 },
	},
	// (2^195)*G * 9
	{
		e: gf.GF255s { 0x39611A8FE1675B2C, 0xFF7E1AF27A122549,
		               0xBF6F1A757B6325F6, 0x074A5FED55883EDC },
		u: gf.GF255s { 0xF2B5C06FCAFEB583, 0x2D0EFBE6A1A2C97B,
		               0x466067CE1442F46D, 0x13349ABFC999DD67 },
		t: gf.GF255s { 0xA4D40FF42D09A54C, 0xC003C77D16DF1800,
		               0xA0761711E6393676, 0x2C1C0AAC64969117 },
	},
	// (2^195)*G * 10
	{
		e: gf.GF255s { 0xA3062A05A4EEFBF4, 0xE5D5C2E0726FAB1E,
		               0x60DC0E92531F05A1, 0x4643A40C89936E48 },
		u: gf.GF255s { 0x1F87C746453A8A8A, 0xDDCF1D005DA27586,
		               0xC77308C7C7664BF3, 0x5774724E6436D94B },
		t: gf.GF255s { 0xAA14EFB76CB73C25, 0xAFD7FD333DC2B252,
		               0x58B8811ACCBF5685, 0x1C9D11E40471A372 },
	},
	// (2^195)*G * 11
	{
		e: gf.GF255s { 0xD13062D63B886BC0, 0x59FB0704B8D4C07B,
		               0x012B93CE6A5FE90E, 0x5884A598A4404810 },
		u: gf.GF255s { 0x08A8CB0C9987EBED, 0x05428BBDC4D7168A,
		               0x5FDD560F0207E6AA, 0x4CD603F72497FC6D },
		t: gf.GF255s { 0x3F08842CBAB59BD5, 0x1502C7440024CECC,
		               0x7241E80F6BE0B6EA, 0x6942504CFB557160 },
	},
	// (2^195)*G * 12
	{
		e: gf.GF255s { 0x94CCC8EE7B7A2079, 0x98E6462E9368F9D1,
		               0xF97BFF023D44A6C4, 0x5396F0454D545564 },
		u: gf.GF255s { 0xB4754DE27F72D537, 0xC0F01E802AC5B1C3,
		               0x6492E0F04B64FC09, 0x14C0F255E1CDCAD1 },
		t: gf.GF255s { 0x392EB3FD516E0A97, 0x4EC3DC5605B71FB1,
		               0x7F75ADBF7BDCED64, 0x6FC86668C96482D6 },
	},
	// (2^195)*G * 13
	{
		e: gf.GF255s { 0xEAD383ABBD6BFE31, 0xBCB0813E9F6E73B8,
		               0xC58E3723BA25DFF8, 0x4F9DA6DF6C6C9E5B },
		u: gf.GF255s { 0xCF158917F580DF1C, 0x0E14A880D3E90C45,
		               0x60D780AFF882B1C7, 0x420C0614ED165A4B },
		t: gf.GF255s { 0xCBDDA904B121A745, 0x31611DAEC9A6880F,
		               0xB43A2ACF31FEC743, 0x5F142474E4E1F3CE },
	},
	// (2^195)*G * 14
	{
		e: gf.GF255s { 0x37DADA6DD884B154, 0x9BD20465051FC1FE,
		               0xF618BF0ED36E5167, 0x6713DDBB85F6BE8F },
		u: gf.GF255s { 0xE3DF626AA182F658, 0x0AEB0A55D410C765,
		               0x645809F39DD03BA1, 0x40C860215AE057DB },
		t: gf.GF255s { 0xB2C0705FCF581D58, 0xB136DC9FFC64B6A1,
		               0xB0115208919623A2, 0x15D34F3283722184 },
	},
	// (2^195)*G * 15
	{
		e: gf.GF255s { 0xA304489B2ECAC269, 0x8482C85F3145A960,
		               0x4A240F132CDBD9D0, 0x2AE8064429DB3D90 },
		u: gf.GF255s { 0xBA4D3B438436B829, 0xFD932027B0906610,
		               0x43EFF53756D1C683, 0x0C078789FF8301C5 },
		t: gf.GF255s { 0x206A563B0F98A082, 0xBFEFB28BC630CCDC,
		               0xE3FA97515BE2309A, 0x71698D07904FD7A5 },
	},
	// (2^195)*G * 16
	{
		e: gf.GF255s { 0x320F03CA494F6306, 0x5AA587118E79B929,
		               0x4BAD3FD43C3A18CF, 0x255DDF7E0DF3849B },
		u: gf.GF255s { 0x2C94910518A991C5, 0xAAC1F71E18F00C14,
		               0x2E1D217E9768A122, 0x50A69209BB0B3E0B },
		t: gf.GF255s { 0xF93D8C083C888F1E, 0x5711285559145CCC,
		               0x1F0F91EC44969A4C, 0x7AB91A184C260615 },
	},
}
