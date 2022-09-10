package jq255e

import (
	"testing"
	"encoding/hex"
	"bytes"
)

func TestJq255eDecode(t *testing.T) {
	for i := 0; i < len(kat_JQ255E_DECODE_OK); i ++ {
		bb, _ := hex.DecodeString(kat_JQ255E_DECODE_OK[i])
		rc := Jq255eCheckPoint(bb)
		var P Jq255ePoint
		rd := P.Decode(bb)
		if i == 0 {
			if rc != 0 || rd != 0 {
				t.Fatalf("Neutral not decoded: %d, %d\n", rc, rd)
			}
			if P.IsNeutral() == 0 {
				t.Fatalf("Neutral not decoded as neutral\n")
			}
		} else {
			if rc != 1 || rd != 1 {
				t.Fatalf("Valid non-neutral not decoded: %d, %d\n", rc, rd)
			}
		}
		e2 := P.Bytes()
		if !bytes.Equal(bb, e2[:]) {
			t.Fatalf("Point not reencoded properly:\nsrc = %s\ndst = %s\n", hex.EncodeToString(bb), hex.EncodeToString(e2[:]))
		}
	}

	bzz := jq255eNeutral.Bytes()

	for i := 0; i < len(kat_JQ255E_DECODE_BAD); i ++ {
		bb, _ := hex.DecodeString(kat_JQ255E_DECODE_BAD[i])
		rc := Jq255eCheckPoint(bb)
		if rc != -1 {
			t.Fatalf("Invalid point reported as decodable (1):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		var P Jq255ePoint
		P.Generator()
		rd := P.Decode(bb)
		if rd != -1 {
			t.Fatalf("Invalid point reported as decodable (2):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		e2 := P.Bytes()
		if !bytes.Equal(e2[:], bzz[:]) {
			t.Fatalf("Invalid point not normalized to neutral:\ndst = %s\n", hex.EncodeToString(e2[:]))
		}
	}
}

func TestJq255eMapBytes(t *testing.T) {
	for i := 0; i < len(kat_JQ255E_POINT_MAP); i += 2 {
		bb1, _ := hex.DecodeString(kat_JQ255E_POINT_MAP[i])
		bb2, _ := hex.DecodeString(kat_JQ255E_POINT_MAP[i + 1])
		var P Jq255ePoint
		P.MapBytes(bb1)
		bb3 := P.Bytes()
		if !bytes.Equal(bb2, bb3[:]) {
			t.Fatalf("Mapping failed:\nexp = %s\ngot = %s\n", hex.EncodeToString(bb2), hex.EncodeToString(bb3[:]))
		}
	}
}

func TestJq255ePointAdd(t *testing.T) {
	for i := 0; i < len(kat_JQ255E_POINT_ADD); i += 6 {
		var P1, P2, P3, P4, P5, P6 Jq255ePoint
		bb1, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i])
		bb2, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i + 1])
		bb3, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i + 2])
		bb4, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i + 3])
		bb5, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i + 4])
		bb6, _ := hex.DecodeString(kat_JQ255E_POINT_ADD[i + 5])
		P1.Decode(bb1)
		P2.Decode(bb2)
		P3.Decode(bb3)
		P4.Decode(bb4)
		P5.Decode(bb5)
		P6.Decode(bb6)

		if P1.Equal(&P2) != 0 || P1.Equal(&P3) != 0 || P1.Equal(&P4) != 0 || P1.Equal(&P5) != 0 || P1.Equal(&P6) != 0 || P2.Equal(&P3) != 0 || P2.Equal(&P4) != 0 || P2.Equal(&P5) != 0 || P2.Equal(&P6) != 0 || P3.Equal(&P4) != 0 || P3.Equal(&P5) != 0 || P3.Equal(&P6) != 0 || P4.Equal(&P5) != 0 || P4.Equal(&P6) != 0 || P5.Equal(&P6) != 0 {
			t.Fatalf("Equal() malfunction (1)\n")
		}
		if !(P1.Equal(&P1) == 1 && P2.Equal(&P2) == 1 && P3.Equal(&P3) == 1 && P4.Equal(&P4) == 1 && P5.Equal(&P5) == 1 && P6.Equal(&P6) == 1) {
			t.Fatalf("Equal() malfunction (2)\n")
		}
		if P1.IsNeutral() != 0 || P2.IsNeutral() != 0 || P3.IsNeutral() != 0 || P4.IsNeutral() != 0 || P5.IsNeutral() != 0 || P6.IsNeutral() != 0 || jq255eNeutral.IsNeutral() == 0 {
			t.Fatalf("IsNeutral() malfunction\n")
		}

		// P3 = P1 + P2
		var Q3 Jq255ePoint
		Q3.Generator()
		Q3.Add(&P1, &P2)
		if Q3.Equal(&P3) == 0 {
			t.Fatalf("Addition failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(P3.Encode(nil)), hex.EncodeToString(Q3.Encode(nil)))
		}
		Q3.Generator()
		Q3.Add(&P2, &P1)
		if Q3.Equal(&P3) == 0 {
			t.Fatalf("Addition failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(P3.Encode(nil)), hex.EncodeToString(Q3.Encode(nil)))
		}

		// P4 = 2*P1
		var Q4 Jq255ePoint
		Q4.Generator()
		Q4.Add(&P1, &P1)
		if Q4.Equal(&P4) == 0 {
			t.Fatalf("Addition failed (3):\nexp = %s\ngot = %s\n", hex.EncodeToString(P4.Encode(nil)), hex.EncodeToString(Q4.Encode(nil)))
		}
		Q4.Generator()
		Q4.Double(&P1)
		if Q4.Equal(&P4) == 0 {
			t.Fatalf("Addition failed (4):\nexp = %s\ngot = %s\n", hex.EncodeToString(P4.Encode(nil)), hex.EncodeToString(Q4.Encode(nil)))
		}

		// P5 = P4 + P2 = P3 + P1
		var Q5 Jq255ePoint
		Q5.Generator()
		Q5.Add(&Q4, &P2)
		if Q5.Equal(&P5) == 0 {
			t.Fatalf("Addition failed (5):\nexp = %s\ngot = %s\n", hex.EncodeToString(P5.Encode(nil)), hex.EncodeToString(Q5.Encode(nil)))
		}
		Q5.Generator()
		Q5.Add(&Q3, &P1)
		if Q5.Equal(&P5) == 0 {
			t.Fatalf("Addition failed (6):\nexp = %s\ngot = %s\n", hex.EncodeToString(P5.Encode(nil)), hex.EncodeToString(Q5.Encode(nil)))
		}

		// P6 = P5 + P2 = P4 + 2*P2
		var Q6 Jq255ePoint
		Q6.Generator()
		Q6.Add(&Q5, &P2)
		if Q6.Equal(&P6) == 0{
			t.Fatalf("Addition failed (7):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}
		Q6.Generator()
		Q6.Double(&P2)
		Q6.Add(&Q6, &Q4)
		if Q6.Equal(&P6) == 0 {
			t.Fatalf("Addition failed (8):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}
		Q6.Generator()
		Q6.Add(&P2, &Q4)
		Q6.Add(&Q6, &P2)
		if Q6.Equal(&P6) == 0 {
			t.Fatalf("Addition failed (9):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}

		// Adding the neutral should not change the point.
		var Q7 Jq255ePoint
		Q7.Generator()
		Q7.Add(&Q6, &jq255eNeutral)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}
		Q7.Generator()
		Q7.Add(&jq255eNeutral, &Q6)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}

		// Testing negation.
		var Q8, Q9 Jq255ePoint
		Q8.Generator()
		Q9.Generator()
		Q8.Neg(&Q6)
		Q9.Add(&Q8, &Q6)
		if Q9.IsNeutral() == 0 {
			t.Fatalf("Addition of negation failed:\ngot = %s\n", hex.EncodeToString(Q9.Encode(nil)))
		}

		// Testing sequences of doublings.
		var Q10, Q11 Jq255ePoint
		for j := 0; j < 10; j ++ {
			Q10.Generator()
			Q10.DoubleX(&Q6, uint(j))
			Q11.Set(&Q6)
			for k := 0; k < j; k ++ {
				Q11.Double(&Q11)
			}
			if Q10.Equal(&Q11) == 0 {
				t.Fatalf("Successive doublings failed (n = %d):\nexp = %s\ngot = %s\n", j, hex.EncodeToString(Q11.Encode(nil)), hex.EncodeToString(Q10.Encode(nil)))
			}
		}
	}
}

func TestJq255ePointMul(t *testing.T) {
	for i := 0; i < len(kat_JQ255E_POINT_MUL); i += 3 {
		var P1, P2, P3 Jq255ePoint
		var n Jq255eScalar
		bb1, _ := hex.DecodeString(kat_JQ255E_POINT_MUL[i])
		bb2, _ := hex.DecodeString(kat_JQ255E_POINT_MUL[i + 1])
		bb3, _ := hex.DecodeString(kat_JQ255E_POINT_MUL[i + 2])
		P1.Decode(bb1)
		n.DecodeReduce(bb2)
		P2.Decode(bb3)
		P3.Mul(&P1, &n)
		if P2.Equal(&P3) == 0 {
			t.Fatalf("Wrong point multiplication result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P3.Encode(nil)))
		}
	}

	var rng prng
	rng.init("test mulgen jq255e")
	for i := 0; i < 1000; i ++ {
		var n Jq255eScalar
		if i == 0 {
			n[0] = 0
			n[1] = 0
			n[2] = 0
			n[3] = 0
		} else {
			rng.mk256((*[4]uint64)(&n))
		}
		var P1, P2 Jq255ePoint
		P1.MulGen(&n)
		P2.Generator().Mul(&P2, &n)
		if P1.Equal(&P2) == 0 {
			t.Fatalf("Wrong point mulgen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P1.Encode(nil)), hex.EncodeToString(P2.Encode(nil)))
		}
	}

	bb, _ := hex.DecodeString(kat_JQ255E_MC_POINT_MUL[0])
	var P Jq255ePoint
	P.Decode(bb)
	for i := 1; i < len(kat_JQ255E_MC_POINT_MUL); i ++ {
		for j := 0; j < 1000; j ++ {
			var n Jq255eScalar
			n.DecodeReduce(bb)
			P.Mul(&P, &n)
			P.Encode(bb[:0])
		}
		str := hex.EncodeToString(bb)
		exp := kat_JQ255E_MC_POINT_MUL[i]
		if str != exp {
			t.Fatalf("Wrong MC mul result:\nexp = %s\ngot = %s\n", exp, str)
		}
	}
}

func TestJq255eMul128AddMulGen(t *testing.T) {
	var rng prng
	rng.init("test Mul128AddMulGen jq255e")
	for i := 0; i < 10; i ++ {
		var n, k, s Jq255eScalar
		var c [2]uint64
		rng.mk256((*[4]uint64)(&n))
		rng.mk256((*[4]uint64)(&k))
		rng.mk256((*[4]uint64)(&s))
		bb := k.Bytes()
		k.DecodeReduce(bb[0:16])
		c[0] = k[0]
		c[1] = k[1]
		var Q Jq255ePoint
		Q.MulGen(&n)

		var P1, P2, T Jq255ePoint
		P1.Mul128AddMulGenVartime(&Q, &c, &s)
		T.Mul(&Q, &k)
		P2.MulGen(&s).Add(&P2, &T)
		if P1.Equal(&P2) != 1 {
			t.Fatalf("Wrong Mul128AddMulGen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P1.Encode(nil)))
		}
	}
}

func BenchmarkMul255e(b *testing.B) {
	var P Jq255ePoint
	bb, _ := hex.DecodeString("40bb85fb77b5bc0729686725ff9a89c749d64471d4e994931e834d6972fb652e")
	P.Decode(bb)
	var s Jq255eScalar
	bb, _ = hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i ++ {
		P.Mul(&P, &s)
	}
}

func BenchmarkMulGen255e(b *testing.B) {
	var P Jq255ePoint
	var s Jq255eScalar
	bb, _ := hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i ++ {
		P.MulGen(&s)
	}
}

func BenchmarkMul128AddMulGen255e(b *testing.B) {
	var rng prng
	rng.init("bench Mul128AddMulGen jq255e")
	var n Jq255eScalar
	rng.mk256((*[4]uint64)(&n))
	var Q Jq255ePoint
	Q.MulGen(&n)
	var s [100]Jq255eScalar
	var c [100][2]uint64
	for i := 0; i < 100; i ++ {
		rng.mk256((*[4]uint64)(&s[i]))
		var k Jq255eScalar
		rng.mk256((*[4]uint64)(&k))
		c[i][0] = k[0]
		c[i][1] = k[1]
	}
	b.ResetTimer()
	j := 0
	for i := 0; i < b.N; i ++ {
		Q.Mul128AddMulGenVartime(&Q, &c[j], &s[j])
		j ++
		if j >= len(s) {
			j = 0
		}
	}
}

// Strings that decode to valid points.
var kat_JQ255E_DECODE_OK = []string {
	"0000000000000000000000000000000000000000000000000000000000000000",
	"8d94d29406f613c642019d7dbda7121510be572d323181af1a0914e8d1080903",
	"23a121249d7a31475a8603bfab1ec71d94fd0dd033c4a6111937ce8649474b51",
	"a040d74b0374447b31b3c70c815c380c5ab9f9efaca7ee75a51a2a7fbfabc615",
	"892180c7ba29f4980fbd56219704b6eeb161053867294aed2cdb520301721a61",
	"8eabe2a42cca04c02a009faab48eb9255c9594b6c9467b03479c00c22bcc423f",
	"ac3bd3c64e2e441b19b7f35e450c6c68346e76bcb2a8b137a7de9e0185bd9603",
	"e53a91660b5e8e8cbdf8603af5b6e4dca4b70d400739c292cc85ff3f4ba2e875",
	"8d33a31606aec0eae122d81896982c8866a1a9f1ebd750d944e83cac3f2f6041",
	"39415ca52c40d1233087280b99a80bf946156cefc5a0fcfc4cd6f640dfbcfc45",
	"4a7bb7961c063e7e1aad732efc6f2a4ac990318b0a9b406ae3c20c8b64bef239",
	"36530914a8e71d6092ef0414e7adfc2485870dd3b2c967b245a0de3f11f84941",
	"36170b60f0fbfda49917de84b99c3f3d901c7fad3e4210fff5f8be7aaf3b1a6c",
	"2a9d3fa655934f2909de0a158ae66e35b268b2682c795e77c5e9d3cdecf76142",
	"e6eb5e9ed983e1f270a72bd0c4a868798ea042efc8bda3adb78d5bc2438b5411",
	"f153b63a831d9a5309dc6936453883c31a929af9943bf7bfb35c10fa9f4ff968",
	"89aeb7a251b63d70ece76dd1b16e14bb97a6c181c880d369e373c4d2854f9315",
	"bd66a62f55d033b7c9d8dc68714787b1baed252c66220e0a5d6960ddc8f0d362",
	"ef782de017160be60b3bfff454b0fe589c34a3e6f1f9275b520a91c4ebc85945",
	"faaf35e1911512270789f220d2d31e21cf4b2ba6b98da2dffb5fb0dbb06ae402",
	"8f4b5186669f788374e3759badbc10384e9b1661090c2183ee19b4990056b15f",
}

// Strings that do not decode to valid points.
var kat_JQ255E_DECODE_BAD = []string {
	// These values cannot be decoded (w is out of range).
	"26b7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
	"e77f719c0bb70774d38b65f54e686da67c5cb40327dc54a35e4eede259c40ef5",
	"1da164d9f10443a44adaa689b304fc620b0084b83a3a980745f0a7bc909736a5",
	"6685163e6eb1ef12173b1292a5ce168f0902df2d0cc5f8d84221e218455e66cc",
	"9be913019bee05796667cdfa16f24cd352d966dd6f7d0c34760e6ab6a72aa4bc",
	"c7d0dc80aaa13be5c83eab4a2690c8689a86b192848e0b23037e6432afccc0bd",
	"ddb2ada3211019b1f9360c1c773a935b3e0dbb554db3e3296b3edbe31468cdbd",
	"f190aba112cc538499ccb40d0a03981a357fcad60bc6a0b76e2b01a3189c09f8",
	"d689e881ac051dd2d853d91a2abebdea4c587696a207e05e77c1b4228bdfabea",
	"4f62d302458e99fa0d7d3f09dd853fd341414db3c0c84bac9c61b18edfc5b8ac",
	"fcd17606aacd5f2d45360a26c8709680766b6a3167815ac591168b16bd299ce2",
	"ec4ecc7c8683f7ff435012a8458650e4edc3faef26dd02c455b05d659db194d2",
	"b85eb4cf83a14cff1c6360d619a61062fa2a12e52be9d03ac39a4baae25b87f2",
	"692ab83e0120189cab03921949edc4181310f984f2d611f53a1f944c62dfbba4",
	"d22892e900c7fb06a8b481348e56cb6c6ac7dee630ea8d7750cd5d7c809108f4",
	"8f275e14c648647cffe8968b8f463a9e6625510c13b90121e8679a4f18af3dea",
	"38cb8c5808478faabc4aa9bdc27469867d968b3db335ad5711971bc7e23107b2",
	"af39fab77d9de6dfe58ab9f5973588211fb12b6213ffa66038b07fb571c85bb5",
	"f16d6faed0a0acf357900acb15d8b27d10dadf5db0b5fd545eeb8b4784e2aeed",
	"f3af66c2a8f69db459cffb0f18136c08751596c03c2380bd8bec413e758ea79c",

	// These values cannot be decoded (w matches no point).
	"464b242d1c4c8860f05c36c21568a428433c32efe8634382bdb38ef392527840",
	"3a7f1237413bafd046210bc6974b04f8e51c944eaf6b0346c9ebde1dd3dfca4d",
	"840b63fa443bf4f8a491bd26abbfe415fe0898c176b3cc74f54494b15d01fd14",
	"38c4ef0dd9b16de9f38f2915fba12ea7579687b3d045680aab2781abbb54e476",
	"b0215102edc6a067471dfaecd71f7339a08ba005461ee652f2eee219806a330a",
	"7ede36ac262d13c1445e84c7cb0643e3dfc81d0d2b45a8efc0c5eb0fef9adc65",
	"225ee2a9f8cbf55ba44529e5e06db708b0bb9f0724654120a7f642d9f90f722a",
	"1648b3502c9406fb650ffd5fd185f7f04c7a37183e024345e4c060c11da4334a",
	"877c919914c70af38631f5a28ee98b2c85c0188facc3158a44a529f840552849",
	"7112294bcc34bfed1f4050df27e669721d4943387893c806c9b436f5e7666711",
	"e219124c9afa4f7c69b2a010c6a9cf82a8aaae67d28d01c274432f9018d2dc78",
	"98a7fa5f994c06c15b8e919c4281459e1d3f9f6477ba2589346b11fa97088e73",
	"73a881c422e9594ac9d8fd3eff56b8a920f50819916c96aab5b0da0912d75f3b",
	"1fcbb61c5c447348aca7d005f605787950a6e58d91ab4481ce88bc93e7771921",
	"886e412db5af6173f59696d486d6fd80c302febeb740bf64ced4dfdc0b7cff7a",
	"8d7d3d424e413e7e0d984b518dca099f0a6d95d7335ecf72b369c0f4fa12524a",
	"d9321c34e92d4faff2e22e6b45883bd1cc02f30c21f9c36924e4ca9e481d994b",
	"7c33da291f0df8180c031975041b05d684ae46ace81023ff1ef36d137e174328",
	"0577366515537e4d4129c358306105508d405cb310888a02b0361815aa91b769",
	"b17a73c541514d26a361f23ffb635e3ec4d671deb26364a8b608c9b4e200b506",
}

// Mapping of bytes to points.
var kat_JQ255E_POINT_MAP = []string {
	// Map-to-curve test vectors for jq255e
	// Each group of two values is: input bytes, mapped point
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"723a35545d0a7a28943ca04c4b617ed4dcbbb7849eea0a44559756bf8d45d680",
	"f3cba9554c4c06272bfe3f51261c7cb3678c527495137618530c74da2a43414d",

	"0e285c3b0dd3bc06de5f84508a82548f95abea361cd851473572343d5f1b2c41",
	"9c87035855c71ba2691011b35e0dca62a2bc707a532ee80dc64c40e80e828e1b",

	"b09e6ff542a126b128d0c959a7dddb74defddba36804d37a8ff618df0090d443",
	"54d75a829df2b19a2a4d91e700ad918bc02bb6f22107b4983bc1e05ff3190d41",

	"6752e27a0ad5286a189940663c8703891d794d225c70a43733fc7121a9888026",
	"1530d881b5d802744c6994b4244b9cf21da2701d6a0a7a8ef0e29b9f4ff63f74",

	"1f32b5d73288138bb2bf2366bd22e454ddfd8c48feecfbc33ef840c04e45af27",
	"a6dfea8c3442c5ba1fc653f196215590b025f081a538ca51d72869e7ad1ef343",

	"7bb73f7b85dbf165fb2c185bea9bab794730ada7ac12787785555e3b628325f2",
	"ac41b778c468b8c939c251a03dec03c2a49cb82d672a800dec8eee5a237b9958",

	"73689d8361ffa3d866db43f4f07cdd3cc40dcfeb4c87a7bcd0707caeb8b87d7c",
	"4e7960431a67ea2498e8b5bcb515c56de6304a1f486144961ac30b242a16cd33",

	"e69f4b22d05dc7d00763063fd29c4dec1d3250f10d245087549260f257752f9a",
	"73f89a33902e496c5d61ca609fcbe0b385b5ef0357420f23809fd60d703af603",

	"fc911bbdadd2ae6010b5fdff3a236eea7270f63977c6f2d7cb7e404131941e17",
	"b40efdd01d6b1678345352d17a93596873cba4f861f7febbc75dc729b0b3a967",

	"5cf1c3c547c6b6c35a6f1e0e18836ba9495a1340eb650795a3c0353441988e49",
	"8319ff7a59791640b8d811ad0903907e4ad092505dc4f8ffa6431cd6928df844",

	"aa147bc7eb3659dfa327f0af067ef8ebf233e8b95633a26f97da175372a81501",
	"9d3e1587230c097098229a28487ffbbb2ac7007ebefe3d2be631b9925c9f2901",

	"512729297f0cd62a02a7ca546ab8666785aa63a51d0e2edb1d1cc5e760793907",
	"e787a112ac61bc5ce399565ae6c37225a424834c87d2a91dfe345a01f688034e",

	"25195837a7a1fe5585a9fd125809fd6caad0217d01cc0707bc2591f0dd406932",
	"756e0e09b3b485ab4c6b876b4be36e3ae2d57ecb9afdf317ade3f2bd21ee3507",

	"5e78ba1b1d0f5c006b5c0e740d8f37b3f3332c51aaa876debd169d629d313d8f",
	"52595e618a146ee56e11d3a2aa12a7746c57dd8bab13f56f122fb5b425cbfe68",

	"abd9fbee5625b2b9527e6fc007fa8e5166fc4c3b107486fedf27ee9a54b6fd60",
	"26a742ac2bc65720790960915af5652952356d5ee84cf61da72f67a064ce967e",

	"c25163787f737f1e016a67d4dbad8cbe4331787cc65f23d7383f93888b2b9d52",
	"3a214343560cee77d9065f6a1f3aecf523e3caec815c39ce83c70441dc03cb3e",

	"39f5603a06d2958dd485c9862aedf0bca0e71dee6fb1509c1fc51bd457fd19ff",
	"7ea9f6b375a5e4877aec88c0a151d6e7f0e17bf890d63a3ee798d346f7df305e",

	"274822b4713f0d03cde0410fc70545ee608aabd2adbd671495e1737b5f4a0ba1",
	"6d2c6af450083f10abb853449e3a02ed946b4eaeab83b24d031c4a66d66c975b",

	"320ada54aff47138c3a85852524536b4f8817226dec94d9b48ad3b0bdf9677ce",
	"352a4595e6144cc3ebdc5e7ba9d2d5a70f52f4f07d0580cca42304436b815642",

	"2b626f86360f5e1b9f8236e3990d5594b477bc52e0465777be80cb8b021bf8b1",
	"e1b3b53403a843269201e043a73916bd0be4577bfd6b78fb941437f42ea4dd7a",

	"e0549613fd58b01c1e32404c83d390c6acbc1cd69fdde89bb0f4074723f734ad",
	"b4a49c5529ae2f3ad3a7796f95a00ead0969243558f3b6ff6a359fc21191490f",

	"6b1d7b15aaf41ac151923d58a4c8c027566cec8220d98cd4501309daa259d605",
	"1451d6ea3f3c93dedc3b3159c33c3c84c2c3766ee1d3ef96955914bbdf15f069",

	"c49b1777775059bbfa05fd83ec6d14d46cbef91fa3e9fe9a2e07885d926947d5",
	"f9d09fae607c55af5a4dbe073f1cc19f034f8ac1106ed9279c0bb766896a5f5f",

	"134c2597128a8867f90dcfd1c33f3ffc57638534d0e3491cc43749afe7be5fed",
	"ab2473703c6f3b602539b9a6616e8e5dfabaaccf39a4a2c021423e1f2b292104",

	"1bf00cc0ca59dedccac2da93b5a3fec9ff65acd0ce92dbe13fb53d6fb688833b",
	"cd452b2bd9cb2128e29b37627865a0becd5cc85f7fd7e5dc61287d33698a7158",

	"0ba02b6a7ce2e2556739299745ea20ffd4df8556b983bb71ddf4d2cbd8aac77e",
	"b934730eaeb0851d5d4db5d2428d485e029103b88c68e6fc9c24af636f86ba23",

	"ef50a39fdb8f88b94f98877cdf97800dbdd3e95f93415b924d985c846e4a5754",
	"e1906955c28b5034e2dbffc4b9429b9540a7441a094fa453bc97cc4fae500a4a",

	"ea5eea9023cd657d1bab53b4a709e961f03672d975154cb7f0b04191c666fd3a",
	"ed99fd5f5344c641326ab44a1670e1ef7af9cbf10ae4d948556fee4ab0e9260a",

	"f1d6ca77d1ad6b53db5f55f97628879e6dbe819f355d97fa3077d899e98cb6ae",
	"1d348e1c2c68c9ecf1494db397547018c4417268e3218a9cdf8bb30e6444fb58",

	"5793d8a97e650da4e9ca4540eb52fa565dab309eb0e9bafb8659f7d47098b007",
	"3b9c9c2623e23be1b99c3169d7b08fb73ca892962233a01cfe1b86e2d0fd9903",

	"0712f766c5f6db6ccd926d92bec2d6a35391ccbced21661acc6ec57d6757d5a0",
	"c49595a2bce3a434d8701ce7bfd8d5241912bb032ac250afe152122e24fc2865",

	"ad80519e5b13da5e05e0652435d944656a28af7734bbd05428cdd8b99609e9d1",
	"fa104d70219c53df9df3dc2370d6a0ace0806068e16ee803e3e2695110396614",

	"c909546acd3ff1d750ae50b886dc83d2e9034a3972d5de3891825ae41b60afef",
	"ccb7ce67b471335891b82a1301d7e98fe8c30d4972f5d7b1b22b418f2eb1253e",

	"abafb6538877ab38cb909307eb11ead4e11ac2054110f174abec506276c30559",
	"b5e61c01ce7f0294e7289cedb19c973ef0ee600285bc3626ec02e40c1293f565",

	"c266be61b0724dad07a2417c78710767baebf168e71984ab084f808528a3fd36",
	"d1a0eb95a8681eb0bf72cc7426e708cbe0acc4d1fac64a0f0b0e70841894df08",

	"ec7ac0c0e996296d3e5c19dc7b14fb3b140cecdfeb7909c8afbc7441b1a25851",
	"36cdb3ed03e9f17378f36ce5b851c86f31b594204e42873329def6991548484c",

	"4b6e7b392d706eeb04c76e552d0ebca92b9155770c169dc40b93fa9b6568c0bc",
	"451bafd5654f097900ad0bc3f860937284187fbfdaeb0ec5a2da083c3d585030",

	"9c31cf29b6fb775ff36abfba52a800ca75ca5ae9ec70e6e334e72af0cd4f6d44",
	"13b10053406de3298cbb83ff7aec300dc66f3b381edf51220a4fe59879daca0a",

	"f1969458ac27c6dd5bbafb01fe0b572b7a2f099737f67289f6124134f2e35c40",
	"08c2bd77a2793ae156c89043a0d41a5dff7c3cc33fee347e7f4a987b0ad0c44e",
}

// Point addition tests. Each group of 6 strings is P1..P6, with:
//   P3 = P1 + P2
//   P4 = 2*P1
//   P5 = P4 + P2 = P3 + P1
//   P6 = P5 + P2 = P4 + 2*P2
var kat_JQ255E_POINT_ADD = []string {
	"dfb1ca41e8528d308ccd8bcba774df736d45c046be7ddbd16cf4d09e4bd6713d",
	"5ca1f4dcdc62691e34a39ef351d89b23279bc0c9df39ce041ea1c1797c2dc954",
	"1965427d058cb8177e90a5347a71c4d8e8c984efb059b391b834fa020276336f",
	"7f95cd111ebcc875c2901d8cbe1a998e515da8dfda172afe7ef831ab7cd9fa34",
	"9dd5baa72fb2a502a92ebce19770953b61edbbe918468c8de917dd0fe5d45333",
	"979158d8691c7905469f7767a093250277fe35814d0d7660a91c701c8ea67c4c",

	"c314c44cfdedb72330565f3484f0b9ca5af5f8bd07c435e295ef236f25d41a5e",
	"21107c12f9ffd2fd3add888ec2db9b1c90dd2af8e3c85d4de5188e711c4a277b",
	"469506f2b9102b9749cfc3b5ee3bd40fe5431671095d816fe4c532a8644c9b52",
	"d9465de29203ad90622c4effd6ff37c99298e7a973742a7ce41c926df2ee5c61",
	"be51fb71259eaa489a24af77e2eeb319dc57b9e5c5fd02fc573eae58e3c93001",
	"9fcd219e21e802f9ca690fed8c88136bad85243b3e58524f7ee0ff0d5a87fb66",

	"1608e067352acfe93e4d9925860e7fe0036cee802a41187f138c7407148b6249",
	"8f7a3119ed80d0f37fee244ce0ad6ee64d9ff33ef96c870a87dbc0e5405a9d46",
	"b0c6ebe2dfb85ac3d7c707d3999728cc8dd835ef8cb96fd139fe6e1828841357",
	"85fa17f27db6337446b9285f037cb4ef13da278189864d5c080b375915ead92e",
	"c86e83cb4da7e786ea6b60d783c8ee438cbeb6d3b48712fe7a05f4ab5799724e",
	"277703b788d7c22c3ef83cb1dadb8f12c2cde24273aa347e93cb53039e0dfc36",

	"9ed7b348a6c89cf4ae7686ca60d7e8e4d84aa8914902de9b0c3c0e773a5d2338",
	"3ae19805b28e8fd11aa08e2b3504a5f929c7d195108e9f1feb6151b50c16b002",
	"b2b6a406654c7605a1ffade8ede158e306e44ebc6c9df121be93d386803d713d",
	"fc9f740da4eb3f789ff9f2de63566c1cfe82b7e9469312f92e5f06add5b42c01",
	"9a19d29014ee1166eb8075802e473fce54ed18ed64dd5d177b3609579314735d",
	"5a16f010530c9e34358a20f5da9a341553fef73b8d0469f025586316c46bf519",

	"e13b4bbc01698f5ad066af91836e421fd222c6ff9155ecbc8343eb596c96625e",
	"304dd928342859ecc0b8f8f8a04653dab139e3b5c19dc6d94285c5426e0ff81b",
	"43b646939e1818e4233641697b95a0a498d14b94c75f2fa40f87ec5f1068822a",
	"b8864ffb961f3902316463cc1c10242eb21783807dab44155e5ac0cf1c79671e",
	"799f6b2a9765582995e8af3c409d662edcee092ef8b90b6fc616753d610be611",
	"e5aeb24c1e7c70a09405918bff754d00514e21e78b30bfb16f4afe0afba04026",

	"456e102ab8cd3de4da2f87aef97de56322cfd090cd1bdfce2f7d1a12f3044c03",
	"dab7d1e828f3a0fa442813138b9500d5a6e8593cc2b258c38f4bc2f6677c2b4f",
	"4ee98f9b18ab3e943dc44ef8a1d34cad80010d56da720ad481d106ead5f8237a",
	"69c3d94898c16ef6a8e8b704bd75f795d5f61a2d0e5ad88ede10f0463511fb4e",
	"adef55541436f51828188f2262d14d0338e29e225619e8b4367373b56e717505",
	"779053e1aa82aab9b4db23a6c62b20c1302ed57e0050f9ccc6c962abba9a2919",

	"47c3a0316cbe7eefd75e9ac28d048ce4a66c4d1f4b6e1ac2f01b92202b0de45e",
	"868e38053880afaedea957f8a15f5fc1a7d259a9e9c81d17c7c6e459a242e343",
	"1ea1e821d1c0455a3bce1b02cba1954222eaec267609f0e020df07a30f4d651e",
	"fd15ff8f1ffb8ab93202477a50bffcb51919f9601f00dc4e8d5f7f67bd2e5961",
	"5dd3d5b70105d7b9489f371a118982ace633004d0d6bd0826b6ce531f1ee973a",
	"6d03d8bf6b7151c4bbcde92060f50cd2079893a4008285ee398fb52ae9ff3f79",

	"ee848f7817d3e971b2b54c8eb7e10b1fbab5ad289923734eb44823a02527f072",
	"a01f4dc1ad54ea5fc196970626b5cf4ef85f70c653b094611d38a5066f2bfb45",
	"11e843d9e9ef7f61e7832326db3ba20c0d32e8536a2d626eb8972288a8415053",
	"cf17d00741fa91e8e3a8e5a9567fc470a37124ad7bb1e3261d42e06c61d94f40",
	"2787fc966f26beb21166505bd8086538d7bab1c478d40afd4aaa6566d1f20c43",
	"6a6a2ed9fbcc30fca854f44c5bbb461f08d8496e51e4448860fbba36131c214f",

	"7a0a00759d1cd4d5e8b105de55df8e24c4c0f77126d93012d821f7329ee8bf4a",
	"896782828d9778495af5e81c1d2728eb65701b5702f5c76339c6327f747f2d79",
	"5fdcd0e79e3710af16a8429971d604a37219df13a81db02f17a075959636aa78",
	"38377262c89d13499fec4289370b84d2f8909b5ff7f78f852fa58de75b501212",
	"33101a47819a4cb010a2fafa2200b9eea178eb295d827aa70ac1eab90f209847",
	"07e8fdc7e3c6f94c7ea00db775c57a78e276a9bdadd9bf48edc8fa2e7d84be7f",

	"20a3145554aa24612918a65010f5721dc35f1b521643bb7a9bd128283eccf801",
	"7936bc60842ba0473107baab33f72323434774cdfd1fdee13d4ec630173bbf07",
	"d093d883ce426874c3eb10d07d742945e54dbadfd3d0987ff2257b2d6e36de41",
	"5d7790edde5459bac7944afee96d11ecf7a09fa908d017a1dc0bfa21f48c6c5b",
	"6aca1b02bfa6d8a3ad6cd8885bae5b67f25988e0482b2d89066ae73240fc340c",
	"3d6f0b6f496aa935d91165d70460f2f98077fb07c8aa440220fcd185a9bbe804",

	"161c59317e19a3563366812fece786e996484b8fd5d1c148ebd7dcdc8e88be5c",
	"7f3715d5ce4db0494b4a3385a04611b18a4d7bd22765a212beba300a33dadb6a",
	"aa27b14f18dcd37ea3cc8d240b643c87c26f6aa76ed7e4fe96557c84e10f820b",
	"a016d56b4104213befa959dc13b34adf4721e848c7e06908c8dc549d873d9d64",
	"a15e3b933bd4ec8f3a81a17249cdc110478cf1f936755f854ae1512df2f3ed45",
	"8e0c587d9ebedf7ae2be0a40670068b7a063ed82d0bd2cf8488c83a9a2578e58",

	"716a5981c982111525789284a3b28513647e27d94a8bb8d070f0a4b7d994844e",
	"ef84ea063af30e199e5d1b68c7d0e2ba3474cafe2b1e68b65a9882b660ced801",
	"a30e0c2d8d09e41927fc91751d9c7b93cd857805fe3ba094cbc0f2d81dd47f20",
	"36954215cc4c68afda107cf709bc02b293951cefeb3cdf255c7057bcd0594f2a",
	"ff26b0308c876b707eb05fdde588efcee433302123afe241a7779622ebb00475",
	"f241b4863d570b3400e052c190a58262d7350d4ea23b8a3a51d46c7fb616d373",

	"eaa1f6043a6fc4aacf323f037b7f2d3fb4e76b4294b9ef70bb917c9cd9254d76",
	"592f8c5cee2995928ecb7a4a5bc4c900c90fc3ce160829c99e7d9bd94944de39",
	"a875ac408e17b91a87f02699cada2e861157f431e3a534ec2596ad7f9b5ec606",
	"1eefe83ba8687ec47cdd7f3f452429dcf7f409471a90591a0448cc1ecf0e5e3e",
	"29a810343e6df60509302e68b2a5d0afd8e0052ce66c96c12f1731aaac9c8a26",
	"c9191c9a53b59af68a329c0ecd9bce33b859f37fa5bf810dd13cdad11b2fd941",

	"030605efbc3bb274c4a57facc65d60def032e2c5fc8fcbffa0132ccfd05fc342",
	"b35a69133de37947d0d4a53c978ec6411be97da282e36596ebd3d7c6ae53af7d",
	"c3438989ccd4966ca03afbf5ae4d9b4fffb63233ff5cf62ff2f89518ab859845",
	"03267f801c066774e3ce7ca37c70e2b011389b46ed48cc2bc506deaf349a4167",
	"06169d3dc96c474d31dedfbe5c7b5ebece757f55da11cd39dffa1d1f5876f71a",
	"f30fed8c79c1c5accff0bf1027b3a22dc4766cb123f86a6e2429fff470d1da7a",

	"ee24e224fdf4cb9e03719ab4e621eef4bce1257ae337c9b55306d78d8801df15",
	"bd80f4653d286e6150247bcc0fc3028857fe1af36d034e410928c0a3e6ec0522",
	"16fa8cdf21208d7cfe1a32b1b8a76fad08e3bd4388559d60e551505b0ff1bc2a",
	"124859494ea846482a104ffedac43faaa90255fa16c4e2ad6d56f2e9d928b40c",
	"6fb267fc018e758ae323c287639c3bd070128f060417cd9120d0f050f8015d57",
	"7ddf2a10f4a37299caeec5ccc9a0d601e3276ee0ec3a280e63cc8e777dead80f",

	"89ce0216734a92f4c20f35baede5961645da4429a90ed58380e662d1efc41213",
	"419a372d235448fbae16df3a8e2cb663372e8d90cab86a4cd6644f864957a22a",
	"fe7b16216610acad06fb1b94f0d9ffc61f8bcb54f7dcb3cf51c3b6ac194e6156",
	"4d49806b7bccc8be2aac9cdd76ebf704bb441ff57e714844ea7eaf0a7e19e77a",
	"fb199c8b717288f1d3d4050ffabd0c074a5a3958064ec96b8909f5d55bb1eb70",
	"d01b9c987ba6049327c96f9b2b5dc6b83cd26044522b676c1949197288827a6d",

	"67dbfa93752bb778433352037a2a9074442cd4d245ae48cf1144606cbb0a6752",
	"79611d1fa2eaee38bcb69e015a1ad7b76486171f1730124c8b1bb2df8bcc5218",
	"4f66724b154d19e0ba806e7877f668c9fd36a5e82e5ea42f2f7d887dc8ea5d7f",
	"eac67bac9ab914eda3279f7e6c1b7d65497899b6254fade7eb06d5e7f9fc9615",
	"411b470ae6914ce475c102cc3a4137b164cc261f7651ee66c1eca8d3aa1eaf56",
	"5300a7cd35274040890994d59131dd845c935db8b9dddb073f6fde47f4633430",

	"fce47a5b4d1663b9eeb8d93012096f12a6509f45a7d2de5269b6b12686492320",
	"b607242257bd8805c6798b11b8ff30b3bfa5a186d0bca466a938d2f01bcbbf42",
	"9339b879a6afbfabc0d3cb8477b639d0d2d5483782b4ab92199677cd15e76647",
	"995456185abe9a56a94027a6708f7931795fa8509454a2d0b5f4c87d2182bb5c",
	"2bfc84e851330e272ec66f518cfe84d7d444072b835c4f4b8b9c55dda8b9dd62",
	"8950242ea5180f4a6c6aaa57539ce4d861be902d81a3e2668a0b62a6352e6b58",

	"9af6d2b370ef3d90d1fedc3dd4b77ecb895416a51ac291f42e281ede5c80aa2c",
	"ea82dc22ee1cab5b512aed06024124c270088f54fa257d3a4ef74e311c7bbb0d",
	"44729e9aafc33a1f1d981b86909395f3d5b2aa4da2d2289e0df06232419dc35c",
	"98837dc918083bb59db02f472aaf329ff678a92d96a4b26106ed353fa4a80248",
	"4a0d207b569516ef57f8d11ab777b1cc75ff0e37aec2439bb43b151215bb8852",
	"a4732968e39260f27f1e2f9e91bb29438ab11e8ae03916c914f64bb83864b128",

	"5f65a4910e32732bce8611479b2c64276b426e935981cac8071222fbbec93f52",
	"d382bd0d89a9d52dbc4aac37ca3cd8dc7b656235561baf43c4d9894ffd498f0d",
	"d88df2e1d8d02ba7957fee18949be935e9b28b4e6ff842a020255a5a8e1d517a",
	"4df39c13b5a53ee6cc628f0f6aa833c4263cb19e94790fc9f9d9cb7e701d781d",
	"6c48a7564f72c814bb855604f64387a88a5cc0489f26399b1e646dc4d93ab642",
	"75d6d4dd5bbcccc2c881ee126812a34712e7be5f7bcb4111751aeaabdfd7b23d",
}

// Point multiplication tests. Each group of 3 strings is: P, n, n*P
// Scalar n is not reduced modulo r.
var kat_JQ255E_POINT_MUL = []string {
	"ec2072da1514edb5a2ca0b7febaad86d8824e8fe1ac08d7d2ae5f76cfc5f7006",
	"25fad24c31b6110d4aa6ee0c5c5144b7f13a92e26aaf0e2b8e8d16a36ccb5a02",
	"b3299c8f4ec364e5c5fa06f7a88d8672a6c7636c369c52a4cf514068f570521a",

	"107b5defc23abc386f2604c7e31fdec74dc2bfefe1ebebc1796a432d62a3cf3e",
	"f0362069ebc5cd5720f0d2c909038869b0aa4e9e7c20f55da62f632a525fc251",
	"db2515bf68db620432ccf7ec62b31651a8128c4ff718a03dbfff68265cd1740e",

	"5c29b647c9b50a42ee34f6e7da5baa983aceb7e047f032ad6e05a756ad2aa36e",
	"8df8ac4a59ce0a8ca690adcd296cab625e6c566ca5661b5ed13359053c87e8dc",
	"9a051ea2366ecddb51d160c64dcb2b1d88717c3bad52368edfec5d97035f4472",

	"9d6af454a16fc4278ebc66bdeb44f79b756c35eb209400f43bfb97ac81c49078",
	"2b215685078883dc4352e18b3c4383515203a7f9e1ce3f88eaba91b7fabe0cf4",
	"61747d88ba8462240c2df91b96fb93b668e10706f8a038a4819cdb8e698a346c",

	"4626e9ed341131b33280e5a119415eaa7a198565b98c0f0da63337fd3db4214b",
	"b419ab0e441a7d63387b4e755d404619613aa78f8d220f87829bb714ae7f6e54",
	"f5c24f70c751b3b66c7f61362ab1a1f64e1d9a6d47bea1836b3d5ba08811f064",

	"c925b2a6c3d13e78817353fc79704168ab2930758f08bf0236960a9cb213a338",
	"008602ab93b248808e011f953cb1fa6a7c695e55b158dc48c44a45a04a94bcd7",
	"ad5143b7da21cf6fafdcf6656596452df7e432e104cdd2f5cb7bb23f0428a35b",

	"d589e8edfe9bbfd0693bdd9acbb122b841c497851fac3b16a392475e1ff0ff1d",
	"29f40a70f4773e15cbad978d7fd4145e283b66c2ade19713f6a40d4dee4c8af7",
	"ba71df9662f7af1af66eab0db0c22d8cf939e54cd86db01bb8d56dcf98f5850e",

	"f1813e5593037e45e43912b3a9ec728a11fca7f2476bea0cf26abc0c64ce3677",
	"03af3892743f62ab17a72ba891155e1324f74226da95be9b1f9ca3ffad2ea99e",
	"e21ba016a4a7e1f347d877849875ac053434cfb49fa0af25bd723b78f7ecf94b",

	"0dd8910a46d29e0fa3c5a41ef681b4b8ab48b69e58776c6bb979d00ff9f9a149",
	"53f4c704cf5129551668b3a03b0c7a460ae4d573c008390ef746ef4e0ca439e2",
	"0a515abf5c8b3f4b0a12a45522a85adaf541179dc1ba5e17725b8d5130881303",

	"a986cf4b1d8fc36dae0726968767f3485acc479679066dba31ea6eaaa1dfb65d",
	"5c0dead41ce62aa0aa0261fcc5f7d83094114b3f1d9dbf21436975a57e485a78",
	"dabda21a821c9e108407fd4329e5c5a7af2cd4b39d65fa45591a24523c669e2b",

	"4183c28c9638ca243743bd110e9b515117e428bfcca60744e9a270bbb883ce2e",
	"f77bb98cfdaafe7690243476f42246b361d952b2eb02cfa4fdef8af21602efed",
	"cfe368fa3d8e193037b3241bb4033ec1364c1bef93e6460fee89cba9fafb8a66",

	"acd8beae6e7b510aeb27fe3bf717188ed28efcb10955aac93c14e9943848c303",
	"116377748e346f146d3e45bf41886c84a0ad9fd0cd2a4f8feffe9a3cb5330eb1",
	"23f585a5ede68378b06aff52063073d2930bc2baf6da63f39a44f64a73cd9755",

	"579e14f03c97bd9a7c441d0f742ca39335daf73be03f9a54cb4fcd34d805cd2a",
	"0bce6231deba2e0b897d790d5be2ea2c97a4ccc3b1ea49e4e73178f156a30356",
	"d2af380d8e8affbdbec3149f1df9985f32f79788a58cbcd2719e612643b3186f",

	"7b2a0a04c3759c227d072116df2eebc38c49c85bf83d4e128e5c0649dcb6ca4e",
	"23779da86ae1ab0c0325ed7804609bb6181e2c810d3e1b8dba3ad9ad74f6db0c",
	"5c35b5ecebb196fe96cb2bdd527e81a9137f0032305742bf73d6a092202e782c",

	"eda214b4660ec5ba5afdce874ada58d135eb2655d1aec476a5a1026ef3f42713",
	"d71556a94cf181e7d8b66c0cf0d110ec30b05e785f94954a3c6a63e5e4edcb7b",
	"b2525ce2f5418eccac97078d516b3c67b50d40d8c842f9a1e5ac73f3f0305f60",

	"3f6d3f73cbfdc880708cff7f20dc822ab52023eb9eb5c9a8f5c21a7f35a0c041",
	"303eb353559c0dc1d04ad64251510471e09fb034fb9d688284fab667e34b58f2",
	"59c1202d39f892e5e0e31de1d137d78ff8dfdf88367b565a892ec74074dbe97c",

	"51c906f04c0b764494467c12b62f67e6d20394a3ee7315ddc018fd69b1dd8e0a",
	"7a4b1618f50aec7225e0cee51a0b5db1a952980bec21a70a01db17019cdf25dd",
	"4022789f4b0e2f1592a2cc3a6ff4a5033d0064df1437f9766bf4604c19c69e23",

	"a92c6cc12ae0f7b46c4f1e62160158b6515f119fda2575dfc92d13dc1ad0020f",
	"b4b1d89e0ac4b296de27192eb7c6cd1321ba0d58aae8fe2de55506cb365258b4",
	"4024d0fbb52db77533835155b24f9ea753c70aedd9bcf591966ce5e369eab918",

	"f47c718ac72b4dbf6edcf9493b1578d5f722e366fccff3ea153308d6f11cea67",
	"1f6bca0e690aa713925f387861767fb7269de2b9fe7c23138cffdeab45dc6244",
	"27267b5f56abd9c1d840a4ae8e5459aa49957d9e06b4235ddb05591b9c8ac070",

	"17013b8bcb7376dd0fa2a42a86f15cf8502fd111578ab86089cda81ceaa94266",
	"da04d1e8f12036124f7f4a5d82a13dce379eaeabf67c8f5cefcb8e4079f4bfcc",
	"2c63bdced2971fc0a6873a1bbe013a6a101f1a0d774e8157920e8ebee208d425",
}

// Monte-Carlo test vectors for jq255e
// Starting point P[0] = (2**120)*G.
// Then, for i = 1...10000 (inclusive): P[i] = int(P[i-1].w)*P[i-1]
// Below are P[0], P[1000], P[2000],... P[10000]
var kat_JQ255E_MC_POINT_MUL = []string {
	"40bb85fb77b5bc0729686725ff9a89c749d64471d4e994931e834d6972fb652e",

	"d3c47ce3a042da4f3da80852a8c5bbbda0dcdf4b1ad51c5da9f746e0dc5e5760",
	"0139e284913ca89a1ba4e4b89cb836a5b4ab757a802f3d010d6980ff00821e69",
	"acb6ef96aadf0feb0f17dcf71dcd3d87a7df4db26beff96d09ee8355b03baf44",
	"7d0bb44498165305f7eaea417a368680e155000c26e6234c0ac9d0f257383749",
	"e13a29f7384827f95964a55714de20819d639069cd76d38611eb5810fc6a441b",
	"70e752db5fd86b017cbbcc94c1b6672732d8402fbe6deece81fd8912ebd9642b",
	"b63b6f33168ea7f9ceb6b52617c0cf5a9913a89cd03559db287c1e1a2d86ed08",
	"ebdf3813932ca6d3e30da8c119ac24dd9e85531c3d7b8c85bbcfe7f20554ac6d",
	"8a47e101632ef7dad747ec052619b62263290e7942e83332f8479ce778d8680f",
	"5217256da97d62b3d8c572fd9cd40ec39c1be3122797775ddaf793ba75106200",
}
