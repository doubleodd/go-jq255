package jq255s

import (
	"testing"
	"encoding/hex"
	"bytes"
)

func TestJq255sDecode(t *testing.T) {
	for i := 0; i < len(kat_JQ255S_DECODE_OK); i ++ {
		bb, _ := hex.DecodeString(kat_JQ255S_DECODE_OK[i])
		rc := Jq255sCheckPoint(bb)
		var P Jq255sPoint
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

	bzz := jq255sNeutral.Bytes()

	for i := 0; i < len(kat_JQ255S_DECODE_BAD); i ++ {
		bb, _ := hex.DecodeString(kat_JQ255S_DECODE_BAD[i])
		rc := Jq255sCheckPoint(bb)
		if rc != -1 {
			t.Fatalf("Invalid point reported as decodable (1):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		var P Jq255sPoint
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

func TestJq255sMapBytes(t *testing.T) {
	for i := 0; i < len(kat_JQ255S_POINT_MAP); i += 2 {
		bb1, _ := hex.DecodeString(kat_JQ255S_POINT_MAP[i])
		bb2, _ := hex.DecodeString(kat_JQ255S_POINT_MAP[i + 1])
		var P Jq255sPoint
		P.MapBytes(bb1)
		bb3 := P.Bytes()
		if !bytes.Equal(bb2, bb3[:]) {
			t.Fatalf("Mapping failed:\nexp = %s\ngot = %s\n", hex.EncodeToString(bb2), hex.EncodeToString(bb3[:]))
		}
	}
}

func TestJq255sPointAdd(t *testing.T) {
	for i := 0; i < len(kat_JQ255S_POINT_ADD); i += 6 {
		var P1, P2, P3, P4, P5, P6 Jq255sPoint
		bb1, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i])
		bb2, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i + 1])
		bb3, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i + 2])
		bb4, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i + 3])
		bb5, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i + 4])
		bb6, _ := hex.DecodeString(kat_JQ255S_POINT_ADD[i + 5])
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
		if P1.IsNeutral() != 0 || P2.IsNeutral() != 0 || P3.IsNeutral() != 0 || P4.IsNeutral() != 0 || P5.IsNeutral() != 0 || P6.IsNeutral() != 0 || jq255sNeutral.IsNeutral() == 0 {
			t.Fatalf("IsNeutral() malfunction\n")
		}

		// P3 = P1 + P2
		var Q3 Jq255sPoint
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
		var Q4 Jq255sPoint
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
		var Q5 Jq255sPoint
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
		var Q6 Jq255sPoint
		Q6.Generator()
		Q6.Add(&Q5, &P2)
		if Q6.Equal(&P6) == 0 {
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
		var Q7 Jq255sPoint
		Q7.Generator()
		Q7.Add(&Q6, &jq255sNeutral)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}
		Q7.Generator()
		Q7.Add(&jq255sNeutral, &Q6)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}

		// Testing negation.
		var Q8, Q9 Jq255sPoint
		Q8.Generator()
		Q9.Generator()
		Q8.Neg(&Q6)
		Q9.Add(&Q8, &Q6)
		if Q9.IsNeutral() == 0 {
			t.Fatalf("Addition of negation failed:\ngot = %s\n", hex.EncodeToString(Q9.Encode(nil)))
		}

		// Testing sequences of doublings.
		var Q10, Q11 Jq255sPoint
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

func TestJq255sPointMul(t *testing.T) {
	for i := 0; i < len(kat_JQ255S_POINT_MUL); i += 3 {
		var P1, P2, P3 Jq255sPoint
		var n Jq255sScalar
		bb1, _ := hex.DecodeString(kat_JQ255S_POINT_MUL[i])
		bb2, _ := hex.DecodeString(kat_JQ255S_POINT_MUL[i + 1])
		bb3, _ := hex.DecodeString(kat_JQ255S_POINT_MUL[i + 2])
		P1.Decode(bb1)
		n.DecodeReduce(bb2)
		P2.Decode(bb3)
		P3.Mul(&P1, &n)
		if P2.Equal(&P3) == 0 {
			t.Fatalf("Wrong point multiplication result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P3.Encode(nil)))
		}
	}

	var rng prng
	rng.init("test mulgen jq255s")
	for i := 0; i < 1000; i ++ {
		var n Jq255sScalar
		if i == 0 {
			n[0] = 0
			n[1] = 0
			n[2] = 0
			n[3] = 0
		} else {
			rng.mk256((*[4]uint64)(&n))
		}
		var P1, P2 Jq255sPoint
		P1.MulGen(&n)
		P2.Generator().Mul(&P2, &n)
		if P1.Equal(&P2) == 0 {
			t.Fatalf("Wrong point mulgen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P1.Encode(nil)), hex.EncodeToString(P2.Encode(nil)))
		}
	}

	bb, _ := hex.DecodeString(kat_JQ255S_MC_POINT_MUL[0])
	var P Jq255sPoint
	P.Decode(bb)
	for i := 1; i < len(kat_JQ255S_MC_POINT_MUL); i ++ {
		for j := 0; j < 1000; j ++ {
			var n Jq255sScalar
			n.DecodeReduce(bb)
			P.Mul(&P, &n)
			P.Encode(bb[:0])
		}
		str := hex.EncodeToString(bb)
		exp := kat_JQ255S_MC_POINT_MUL[i]
		if str != exp {
			t.Fatalf("Wrong MC mul result:\nexp = %s\ngot = %s\n", exp, str)
		}
	}
}

func TestJq255sMul128AddMulGen(t *testing.T) {
	var rng prng
	rng.init("test Mul128AddMulGen jq255s")
	for i := 0; i < 10; i ++ {
		var n, k, s Jq255sScalar
		var c [2]uint64
		rng.mk256((*[4]uint64)(&n))
		rng.mk256((*[4]uint64)(&k))
		rng.mk256((*[4]uint64)(&s))
		bb := k.Bytes()
		k.DecodeReduce(bb[0:16])
		c[0] = k[0]
		c[1] = k[1]
		var Q Jq255sPoint
		Q.MulGen(&n)

		var P1, P2, T Jq255sPoint
		P1.Mul128AddMulGenVartime(&Q, &c, &s)
		T.Mul(&Q, &k)
		P2.MulGen(&s).Add(&P2, &T)
		if P1.Equal(&P2) != 1 {
			t.Fatalf("Wrong Mul128AddMulGen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P1.Encode(nil)))
		}
	}
}

func BenchmarkMul255s(b *testing.B) {
	var P Jq255sPoint
	bb, _ := hex.DecodeString("f0f9c911813c90d671dd55050b359c59daa0d39aef85a70196e99c1c1d2bf749")
	P.Decode(bb)
	var s Jq255sScalar
	bb, _ = hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i ++ {
		P.Mul(&P, &s)
	}
}

func BenchmarkMulGen255s(b *testing.B) {
	var P Jq255sPoint
	var s Jq255sScalar
	bb, _ := hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i ++ {
		P.MulGen(&s)
	}
}

func BenchmarkMul128AddMulGen255s(b *testing.B) {
	var rng prng
	rng.init("bench Mul128AddMulGen jq255s")
	var n Jq255sScalar
	rng.mk256((*[4]uint64)(&n))
	var Q Jq255sPoint
	Q.MulGen(&n)
	var s [100]Jq255sScalar
	var c [100][2]uint64
	for i := 0; i < 100; i ++ {
		rng.mk256((*[4]uint64)(&s[i]))
		var k Jq255sScalar
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
var kat_JQ255S_DECODE_OK = []string {
	"0000000000000000000000000000000000000000000000000000000000000000",
	"b1c26ad124f0b9b235e811fb02576bdeed7f6f2bdb9625ce867bbc2bbd751e59",
	"aef7d911554e0e1527f49d49d303e7da00532359589f1d56528005cbac70b153",
	"318b64b9ae83047670e6256054e6b5c0ffc3f9cb24659d809d946a7cf72a7521",
	"66a968dfe975050781f837e3ee6bf0352a969e53e4f66d1276b2c433568ea26b",
	"285611df1d4459f29fe7cbcbb98c752306d75664851836e50c347401a664f270",
	"851a83fabfc063883c88904475b53dce04d664f31e0fa678e297e3d2dee58732",
	"50ce42f8d53d0e1b79185c5be821578dd2a86b20084a6c0eff47b6431457ec64",
	"77859488c089928c355581d2a8e86154b94373ebcaaba991ad4ef19e4fb0352a",
	"bd9f4067b33ea224f281e8fb8b89c9317151796b415d7e491fb469c7289ed058",
	"3a89543388ce55e915cc870a4db2dd4fd9ad2e618f5a5013a45b67505859957e",
	"1bb1041d4123f462baf68ad3e6fc2aee7b2d30291ffd478379ead57470b1281a",
	"8d902970535a38b5e8d283a55d8cfe4ad62291916b000e561ec78e81d30f3160",
	"5758687972e24623a80af50fb524ecced34ade656a701644c764e30517a37713",
	"7f013fcabb9d8d29c2d44764703a4d69b1841db608bb1043e51d39dbd29f8739",
	"9afcc2b95049ee0feb5787056d301115e0b4d22b8ec82cdc4c37e80a7990b901",
	"d696e2470eb8a2a9c86a87791ee99852e4e51cb2be78e826e73fd73037df3a00",
	"510467bcd2af8b15e4833126e9ab0057c86e743420bf6681a664447f115d0118",
	"14ed9c64000ce1100827fa3ec256b7b09d0f3e009dbaccc17ad4e08d62172b37",
	"51fe4489d8c6f777feb613dd9cfc5160a12cc30c2e48bb1594e32bdbfa8a5f49",
	"dcc3f03e68c699d5be18958e90d9efad3ffec41de365807633592ff89d306771",
}

// Strings that do not decode to valid points.
var kat_JQ255S_DECODE_BAD = []string {
	// These values cannot be decoded (w is out of range).
	"8ef0ffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
	"f40f8584aa9e36f1be98b045dfe8573213557185500de9519c317926c9fc48ed",
	"28e123fc26680524a94a51d3e9a120a39e329c383c41c3e7ac87f173c722f0ee",
	"3319b08127d291d4480b441b7489886483ec88eeb77c165bd3698c0b305c0dc3",
	"1895a6aad79f9738f2f6e925843737a9ffb01fc803f948a217255e35cf91e3a3",
	"300ed3309d2df52ae6773cf378be94c961cd1cdf34480a3604efd3840f56f9a9",
	"08572a152d45763d1def57969530cd09b9a52ac7da82b241c0364d081b526af9",
	"bd3bd83424f7ea3f6cba33e25ba93ba0c4c10b6645508048c9e6e467ce1578bc",
	"31d9540603feb90ffc81e6f94975b077ab88448245f78392d52c56334e998aae",
	"7064abf64fc84a0667abd7396c30b63406b677ae078df897b7193208d845a5f2",
	"f35b45b66084fe0b497696f517e3ca9eabcc8d152dcb961d6795fa6ba8c0ada7",
	"0718ea4ce976ab9f64f6c62a5a73c1406aeb4c3c7be2d6b277cf07a8f3b1e58e",
	"25db3d20f89dd4b9cf111912fc3427092464292c8634519e73afafd3821e73d0",
	"e7b7af235613713fbdcfc825c2d877fd3f87df7e83191134e249906060e36fb7",
	"458737cec3f41fe9506af823fadcf0271dd81fd656038f0c3dbbf7e93b98bed1",
	"98061b024220305afb556de3a2becc081d40cdb2e32e0f932f61dde8d668d8b2",
	"1e2f702bd67d6ad8710fea61dd8f5b4e96128ad528ec66f697d73e54d2cb46d8",
	"2da6629388de5b7c5a6ce5a739f8e2f0e6bcec00465f083412161dbb871cf3e9",
	"584eaf85c8017a9c8b40974974f423d882f946ca91f9dc317a62ca65b6380392",
	"6f7b1a58344d99a20c77923f7eb0b09bee724ece288c9624d2343c962c710fe5",

	// These values cannot be decoded (w matches no point).
	"029db5a37622f1b2e4bedb93f641f98b4b26a06936e80a8beb5a0a81c350bb6c",
	"eb23b93c99b8d786f9772371d35b16a444bd22151bf34f97cee4f11f2b54b53d",
	"2c57cc483998414a721a4736f89312ebb2623cd7c7b04318ad2755950002797d",
	"56196701642782ddccf9e812ea8bfe7cbf430519e1f5aea99fa23ba784d0d653",
	"2d144fafa94c2ec516e71aed7967a3aff354bfec9002751a0a07042f4af14943",
	"cc8a98aeb704ac2e884ebba9f9832d686a8e0237efcef408d4deeff7bd248d64",
	"f3c76f1a6c84603b2327e6a09e7c311d3dff636edfd63a6803b6a18543717f4b",
	"6afe78e4600c57add45f3f763536f2828e1ea1a37de45382991fe3ea61957e17",
	"15ff8793e536af65cfa1f3ec63925c2a504d155cbf857023a3905f88af5f5e09",
	"e99847c33b43dceb1605b4ee22cc8e7f1936cac38840b2f851b538f286ec8609",
	"1b57817ef3887174b39bb3fad25cf2a72f9c9bd948c60efdc5e7953973ffd860",
	"bbc110ba28ee2fc04021e12362889b8762b8764c424df62db2d2a57da351c27f",
	"4e2b59d3ea5ea3140b9e0316f793b94bf6af4470f03dc6bed6a2c0abced9476e",
	"4bdeedcc7050b0ff28c7092af6ef78b5bd0e369dea29453850cbdda019155d7f",
	"522db2987ffdead41a2d731fdb7a09122b6ad9909b1fcfdbfdbb5c3b641edf11",
	"f79ca01c0d421023824d0f081aa5b5004d7222e8b1583ff9697632aaffc4ce61",
	"f1c9835aa95a7783884c13bc3963fa29afea36aaa1114b26184319e05190df3d",
	"84f489ca4d49d8b4db13b0ee9cce45e10a4c537fda03941ba8bce5f6a8ef5a0a",
	"9f99165adf8d45e538b4d180747793473eee05e5a42d5db44af76e29b489c058",
	"74f084660923f4de26538adef1b8424de0c2ad4e942b5e1b7c8af056b939e57c",
}

// Mapping of bytes to points.
var kat_JQ255S_POINT_MAP = []string {
	// Map-to-curve test vectors for jq255s
	// Each group of two values is: input bytes, mapped point
	"0100000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"4f347267b8e1655a0a25cc5d302116baba0ceecc838986132a0dfa5bc74e2993",
	"e18ffe262e994600c1471c741ad8e7ae11abb640c398c8b0306cf6d76add594d",

	"97f335d4eada682f09239a3b71a66c11a531064b09328ad8c4ca436a05e2d816",
	"a8eb52e69e2c473318a61c66bc6303666ff2c34ebdace7bead72ba693387c54f",

	"a07e536399349201a5765a9edc51189366797999dff2181b14a9af0da320fad4",
	"eb1e9cd186459e27cb83e6b847ec3153fe3c6c8edcc21f5f4186a060af381149",

	"8350bf235e37b4acfae66083b52061f53b2894a7db9a5f13881670fc48b4309e",
	"31e9a5352c29b0a8c57582cafe6dd7aa2b094f98c1783849c58471551b51e103",

	"3894717e7ed532f279d486c0dd72dc2f319d86cf0407f7975e5ef0c8f6f2d68b",
	"2dd1f1b53821b86fe6784a27fde7cdba5201a09718659e621126cbb757073409",

	"98fea1a8a6c415f0072281f8e11de091e9fc8dee90e8a19c57af2c3fbed4df8b",
	"f93f5d3c0b7a3513982fd68c95ff0ac27e3c6f41e30d49a06de94e33e0c30372",

	"24bc0e9df5b0f627a27eef1ed2e45115f4e04e8865cb0d3d15e9abd6f8d1be2f",
	"7b9526bdb6f92ed5bae05b00000d02b8f348b7be60ebd44cb251ba6e8a8ba13b",

	"b3e07d3bbdc330ffd662d42c37e4c8274bc8b7704eee118e66f6eadba50d4fca",
	"9ec145cf7c074984124bb5cc62998bdc9fe5ec9c34881ecd6041620b564af53f",

	"5ec437395a7df0f0e849985fabc22ca914c1afe4b974de4fe651169a05174d1f",
	"eb0200bb95e62ccac6df2d756051e909a1831c84707c119c1063501c15ba8c74",

	"b79f138c5eac38b5f07bc0245e9aed2179c4d6f90766956b922747b224ef9cd5",
	"ee48b384cf94a9a920692ec527a9116e1de7caacd75c788866d494ffc5f53f5b",

	"cbb53297a9c53c82d57ce8ba26812c788932d19280b45752cad47bf5b81f7dfc",
	"b04ea34f214377afe8a9b75dd252c8a0538e8139aeca59a80490702ba865045d",

	"16161f740ea4c4a639151ef923d5351fe96ab11fb3e84fc708dffedb4efb33b4",
	"7e7988e6dfe4994546f31912fc7752138a144857052f91c1c5ddd7c1cfd76153",

	"e132790de50058fd94aa558c9e9bbdcb87c0b62feea19a83696a90e8b52308c3",
	"37031786f7a6fe73394f8400027901f7534c2a83cf06230198c6f6e0ccc5f374",

	"6a26aeeb3cdfb7bc6bef9f0a76300ef2609f126a1d08128b45da19cfacd7f497",
	"9e008fa86aa90db2a473edd99fd3bf361b232d698a6accf583c922c2fec1270d",

	"c754e692b9c948d513694919e6f9d417cf42a5565d729cfbfe1cd080f4623de9",
	"a7748228ecf7a3afc6a915656719254483d26a753bafbe73b0a02ee84da1e924",

	"3e688cc3e0875232c152b653dbbcc1a22739287a315909f58e04e55e7d2121ac",
	"b5465cdfb8a776430845e856c5de93ad25f98ea78860a0d9a52087942e607762",

	"e054d6b59905617d42b9d31cfb740b9219701f7a902ec91c0c115bf9b5a28583",
	"2aecfc00575cf8e3aa6b981fd2b53fe1846ac17455831fe42c096218be57c504",

	"43df64360cbbd06bf7247495563fc880887451e0c5701a25bc9284fa2ebcd856",
	"5bef1bf644e7c43b6eb46079bfee57fa9efa1ab9f29eb53be6d627625799fc32",

	"fd2e1e57e40c1ef60a21931bdf2329b04edba8d72bfed9fdc186ef76b5075859",
	"4755c542dcd87ac4fb6c4bef8d3b7ff711d58c18d61c04dac44aa391fe60510f",

	"3cee890ab9ad19dd018f37f9162f78f8bed98c0594e50275192d25df071b8ca9",
	"9e0e84fd5bc9e752b7f205962e177068a7458e66764ac19e7c561001bc3fca2c",

	"fe7fe03fe62bd0fe11aa1b3863f3ad86c04d1ddc391bd818528676b84d740c6a",
	"090624bcfe3aae535d039daa55616a1c03f493e7372d0e1870ef4c3142d0132a",

	"82afef12e6bc1ce76f139b581c3e2bde4ca09fe15280ffe32ca2cce87b5475f0",
	"2bc0e5e125ea6754c4df0f49e070094a79f675f3f42b9f880cb894d1aeb69437",

	"97d2c5245f76665faf903a023ac69d02b0399f3c7c5a325300569a5a605fdfdb",
	"51a503689ed77ffd4ee5fbf749ea59a8305bdacc284d5b99853a40544dba5f68",

	"1955f6b258d7283ac4b10ff2913586a834dbfecce17e5cd491054152ab0c88e6",
	"ab718403165ef2bb3d3b2d2450917312b541689585c0a030b5e3b24e58884424",

	"851edef2f6eabe94467279024d96d08314566c64e2e945d95db73607c3a18fef",
	"0055a3b2e5df1e950cd0a85ce97b454546bf35bb422ddb452ed8beacd664d761",

	"ce58e5fb30564c468d34d3b68f61249a27c8fd1e0ebf51c47f68a2465d24bfcf",
	"5f56b38fd1681ffe0d64cb2ac017c0efe8fe18bb2cbb15ff58c03e05d3066d45",

	"04c0845caf76658619d20cfc73a1b52d559325053ee10d867a2739697292d60a",
	"a0ca4ffef9e12b33dd14bc8f58d809ad47326e1f1348c9561a34e3e7f91e4761",

	"19bd05ed8a0deaca33d5782d223aafa6da90edba5b2443b38c9268f130aeddfa",
	"94da1ccd869ff797701897535c00f7d10b488de6ae082902439faa3af4444e52",

	"3b1202c82b34a81a4c996e337492ff99c52284d4c1b00d2fceaffe4c8e95e9c7",
	"a8911426b9a05cae2bc5c6f0c9d9604408efd82972e33417b76d5ae65d61c21f",

	"a98120c5a3b4f86a77b7e6bcd8f35ae136a47c7a9feb9be0907d043bba80ac34",
	"e18cf3de24ea143ae56c6bfc539c49ec3553cde185f47a1a530322527a9d3a32",

	"bf1351db9b57bbb6ee712a2bbdcc6e2df2d54ad62d51a50bdbb13b2f7ac06c3c",
	"2cf4c9072a8df20da57b50576997fd5521d5001bd390e42bac4bccd1a895da5d",

	"dd4e3e5b80c45573e31ae8e8fe82c9cbf7b0dc1f1139718d34807530aa8447a5",
	"fd6ce57786494f651237004effda2ab2adccd861d8823381de97ed427070da7a",

	"433f4deb8e26cdffc6224b3cb4aa3bcb8aa2a6a08f88a7a40f1b2b386560a3ae",
	"66e801a030461a44e160b204a4c45a1bd61f172b493e42eed35f333617a39575",

	"464ed28003f7b8413cdeb77d8f32e45975b349fcbd7ee4e5e5653065f1a460f3",
	"37599895c6a0b2771cf741010790c7022bf181d2367d38d483ccaf088f80146b",

	"e627cdf1507bd03e5fe45c0f29d12ba71076af1f42013d5e24abc0e285b2042f",
	"6bd6607e348174b5c10003e57fbead121468fc847fe358fa2602c9df6d51f253",

	"7a838bd079ef3e37426ff3bcdac05c3d1e259fff639876753012edad3d64977c",
	"f9f4e18b5d503306164434a655546f02f7e4362e833362bf75b64cf4516ae62d",

	"0a89fd51a9abba6f86f2db7c1cbf4346b6b46b6966667d4542a177bac8177370",
	"15bc61c129ffde2c73b1c347274b828955ca01a1df27770404f44f0b020f7010",

	"26ed1f06e41118d43c874bb29d73485734effaf019e23fa8f48d80982318df02",
	"7e23de36a843af132ad5751bfbf21825faca306207f57269db95f53e489aaa56",

	"6d9175799c8ffadeca0ead71dad0557836905b9bc0d94d332f22e4df25ea15da",
	"304c8ab25fe8c3b1e3c49fdc75e769b48f3def51d4505d280e69df7eeb4eb809",
}

// Point addition tests. Each group of 6 strings is P1..P6, with:
//   P3 = P1 + P2
//   P4 = 2*P1
//   P5 = P4 + P2 = P3 + P1
//   P6 = P5 + P2 = P4 + 2*P2
var kat_JQ255S_POINT_ADD = []string {
	"381c4e1998ea738fcee1d18779bebb36271ef55a5bcc3a5595717576e2ce704d",
	"db479deeb0d8e6096b33dd9dd4e18d385308cb553faa89ce8681df208860a110",
	"e0e467a709f40a935a8b1c5e78a745924a0eecc58576c99866b20e435fd7a821",
	"e02e26cc7657170f7e3ca5e9fd7fb11712928b55fb007df9d62f805c49be5673",
	"d60aacf148e642d5f31a4307300ffd46f6235dee653e7ea611484586669aee01",
	"03955a234b2760184666f43a9f3aa1b94330a37d4a2722cfc2dcf436fb67ad52",

	"60f36f052234d98450e0131f5cea458af31b5e36df99c7ac77e8febfd6bd4562",
	"eec699499ef822479c929686a6fcfb68bd7c03004af40a6e7120e61bf946ff0f",
	"011e0b7c6a3332b8d81a44cb34a4a6e5428d1a1e2f0b68973f65331a2aa37421",
	"075b9051e80d2d25c8f289331e5abda12a9046c8269994fbbcd0408f8196805d",
	"a5d16235f08f5db6c03138c280e7b79fdf3abb1f5c9cbbaebbac00d7b5e58d15",
	"90a2dd8393d8f80f07b323a6a3f71aca392ee9355a285d3a924c0a08fdcdb12f",

	"5f999a844a1142c521b30d7962c5d1cc51f44b73310f2859fba957bbde789250",
	"0e463bcd3d7a94dbe0c1db26bb8a1ceaa02f3bb5a158690959aae534a99d5753",
	"8532065bd0f3442f01b315eb0088f42656959eb082e6ff5064e7739383a1e50a",
	"44d0ee20649420c8514f2b152e8551334df9755f18eed00225c0ea4aaff52262",
	"5e026732df7af8d6265a9021cbccd6b4653353171fe228ab3c3a26911a72b421",
	"b0e8da551c65613e0f8e90bbb987941b55eb1818f9c4747cfeba1f6bd715e139",

	"3ea95fbb3baca80e0340315c688ec8b898b17a491f452e2a2e7d7b7d948d2171",
	"2cc136fb22b80d783f6a8162b1292fc9192cd5be4d719c52cff60e6cdd845119",
	"dda4bf35cd7d5c05e49b37640df81a385e4cd9e540e6a84782a57f53543e3b2d",
	"c11fe3be926b10336cbb8f191e31459c6fa492629c773f5d7214536cd299a367",
	"3ba12da38f381c060ac30c20adb436cb816b01fcf72ab777cbae7525767ec044",
	"61987cfeb7e4c45f742d1d6efd5d04b7d7a10246e591d947791d64da569c2d6e",

	"6c44d7a7569c28f80b12273d5ad2308f1241eb1e42d0b92a67213b7e908ad14d",
	"e9b5567f03125c2357ad2190a852c3dec98dfd919721f385066b28f797c64e6c",
	"acedfb558a43ef11ece45956d15a4ce8e206f101311ef7ae51d89fdf1d593304",
	"5974f71abe178125b92b559ee5025e004c37b29e570e1b05702bac49160b221e",
	"9814abe67e6f98b67c7abb8f192d9177b5b707f3e0eec52a08442f3bb6e80464",
	"9e6725094371bfe29b6e114b29bf12eff0167b93b98b611734cdd3e196a0f52f",

	"bee7f22a49b9e5404162552437d11fff2b773b4111ada2ec8dc0faab67243261",
	"c7c2c162ba7a3c0b8ed64176b959207e396b43b9eb5429c418528bfc1370331a",
	"f0615e44ab681f5f6823b3b441a39fb5f518ba448eacce1f104008771e33fb29",
	"95e9fd29ea0c0bacabc1d4158503890784c3c37f3b35c70c9ae45d94caf58329",
	"02c5f804b8607d5525e52810a3d000cd489fe7b79109239b26af3d9d4085d658",
	"95faa802ad86ba6929fc3b01730daa232076aef3b063ebb6b309cfea61f9e741",

	"19ea6d24f1511a233c198ed0895d3fc76161a3e48f310b15cc8e6730cc2c1026",
	"5755e6da5c23c5f0c67b9b1c7f7e8b02aea0ac07498d206b5efa045489e1d737",
	"096ab36803970aa88bf6000f889267543299016c63ea2c629df4e9ed985d862d",
	"8ad8630a4117f64a7d78deded3400ea6adf8f92f98ddfc85c1b1b5819b5d4046",
	"8cbc2d5148c3c78acca57399073370089f535da059a93415788a1528469eea4a",
	"6bbb60d1fc041b826313dc3c0a5116e71767e252cee93c8dc512516979b4ea12",

	"88b303961b1ada474b2294750bd00f1c935f0e6f1de20b07c4e971fb67527f51",
	"823978efdd5b7b305262287900a4eb96cf5e61c351a13a9e3bd61efa687a5866",
	"dd5c8e2ede2ad115a1edcc72c2e3c264e300fce3a4e3d0f09f504e023ac9e31e",
	"20049f4a29799b829533ffe67cbf70279b88bb19a6ffe62b670914b5bafae025",
	"1b7ade87df630960d5295480937545974949714c66184ad3770542f7377e2301",
	"ea9fd3440c0a193d56e93589541d55f25bafcda89d66986689ee171823951666",

	"add5c3e78cb10a1f9ba69997804838250506cc7678649b483e6299e045497371",
	"7483ee768a54e1f3b19422c9bdab65ed8c0688483cf8fc3f80a08d07bcac0622",
	"c9b187ab32b93e9acbd2600d79985eb6e1b0daaea7dee8a372955bf36beba27b",
	"8ae1eecc6db1a3866c2ace308989a0f8ddbfb9b182da844aadc05148aca0e25b",
	"4c8e8c39fe3370d9fa1797dacc26aa93463bbb2cf8309140c43bcb0e16842e39",
	"92082e52bdd0091846e0be899144cc553379b723f609947ff95c6308b3e55a25",

	"00219fcb3c4a5cfb84c0cfb2e0be0c2fc0dc50be01bd15c0849cd2b4e86a8753",
	"61cba11b9c106484a8d3b677fd42ea9aed3c4e3211cf87845db59d085406385b",
	"5c32d46832a9e72ba1b334176d5afb4b499915479559c99ac900abde03d3a503",
	"599f2e7e822c06eeeb98477975b18664556c9de5070e296bb547a5d6b064d153",
	"4b9d81c62c6b1ed1159db5a6b4af8161b7ec864740bc2bb3e49627997043bb0a",
	"4cd2ca03930823ea167e6cdc738cf4bba68960999d26ef6b748ba8ac1fb9f61c",

	"f881d61f1d74ec291e82eb12b274789ceb22541a1d62f58f26c946aa8262954f",
	"a26b9edd4f33a8d27ec17b0834ab1631198bc7a86d34f22eb88b5b59dee40a5b",
	"9bc3688ede79d6d6fcb5f29576e8f4fe00a47ab49913b3ee3a66ec5c901e4137",
	"716f532220d997612022d61df0f2aa890dc087039ddbe6fa315d89e58ffadd72",
	"1c209ccdc7fc984d4d03b763cdfdf0e5ce81d22eeb710057e80b49f2fe7bf74d",
	"8067c4ac59dc0f95470c06efa63820c3ed714146167216eff0ba61fd0ab12c04",

	"6d06d460392a61e9011d18622610ad4ae2b40352a787c440e458bd67e81bac6c",
	"acbe4ea45fb6c601d0b0da1346ddc2315e851def7555c90ef460f9088cd7114b",
	"ed3af00b64e16034cdcd3b5344ca4bcb686f61af5894f685f353b06140b8ef75",
	"d4cff64435af65e2d73db04984c9b1724614827ad272ebb9fc76646885bf5328",
	"2309d522f92b184f7c290a71823f79b9f0c72020c5b6554b4647f813da77ff3d",
	"2f910bd99152d1b56542fd3965ede38ba9647516dddcc1a5bc9db0aa645b8103",

	"94fed56d5124d271fd3b63cd17f4f3bd6f3fb9c9b3a2f2dee7abacfada3ab969",
	"25debeac62ccaf5c8d5019580c7ef97e78c769ca6514e6873e9622cb74b18171",
	"7413d872bb58181d56d0e288d2993b2a552845292632b45bceb2fd8f3a4ef614",
	"fd3526dc4a70604ac9986f916129cfb5995b8aac2093dd06c338cfebcd23062d",
	"8ded66b6e0bfa56a27a5fa4e763545064849e35ee855d65f73546832f9dc0d54",
	"835d25ab5b57ed33402ab7c6701755117619948a558878193fecee12cd6bb150",

	"780fc3cd9224249d1df0ff142f552629846bb1f1beba2621bf79f85a9d314d56",
	"3175dac57c3a693980287572a2c14d0a59a4f078af8d00a79b38ec982b25730d",
	"da270fc356d1a119132b75504682dd56ccbac891d100e8c794d0833285c1c34f",
	"abe881a298959bac0e8ab0b111348a31c3c3c9dada2ccae15a853685f9ac2b46",
	"bb1344014186f4219716a740867b0e7dfc3297de7b35240ec0cfcccb383be57f",
	"9e0bbb128e78780dfa725cd8934e025097be87672f24b5733b41f5f23a5a3549",

	"eed74e9cbd63ddb8f3ca345ebab1b0ba9cc6bd53ec6096fc059354319957fa74",
	"f921125ccf44d17ce2d8d0f29a9f59e0f73d8a759ac8f7ef9250f4d5fb0d7426",
	"1f36c8e821941eb90d24c6341e75d94aeeb2527863a6a634fca5d38f42663430",
	"7e82ad365a5f308390a0ebfd73952cf924050f03620a1e0fec69797a2179235d",
	"aca89dc458bd0e2ca0f73c09e495b2241c085011ab786f707d8c3c1c2f7ebc4a",
	"5c52e4759076e2558be0b5ade6e1bd53c99514c656d61ae4dadbcd63b0934a30",

	"97dc9ac2154cc71d4bb42fb46afb064924e538a88db86d52397cebe124b4b075",
	"3ccc7026466bae54512a3f4ff7f85e703872b74f517a4b607dd5be5a8412a16d",
	"2ddfb977bb7fdc7880efed6953ee1604de69df616964f8a127437e0e98dba361",
	"7185bd0d698a74852ac7639fc9b955bc722f30a019170f01644a353c24e8cd2f",
	"47229988a81941e122554b34bf81e2a75bbacaf44bf263c716cb2725665aa37b",
	"414b57d106ac19110db4815c8e59be5c47af99b139e8ccf13025c39be320195e",

	"fca688816addd850c8930696a49a1e6fb9e14255b03aa86dddd049f6cce8ee67",
	"53f9d223569e492d78fd4b6665ae8994c6a1f852a6c38a60072a87240062d808",
	"706ba0b5f5ce767b83fa2c539e4daaa750e5bd43c95acd5a1eac3d4afd6d6508",
	"9a26cbf9ea28d031068bcf81b642aa9044f0fda980766d22a6b2c8e479b30840",
	"18a65ad2f85cd3c3f2fdb6fd4a4d88e523ec93b19894ffce6f3fb931740d1607",
	"36522522608857cd94e8dd9814bed64f12bf8a6845098b284457376544afd820",

	"293907dbcfab4d4dab6b9edc338997da0d107f3cb47f3531e75c04d23171ee6f",
	"e35d3b0cea71bdfa31be71d19eb58c7a86a4566169f17980d2b32a95b9f3ee54",
	"2f78f90f07e8d242bed8fd5a160fdeda7924d045b77b228ae7d2b9a46e1b1e07",
	"68564e3f97a87ee54e66bdd31f98efd95f90c0081bb348dd179ed547d4732c55",
	"5005f9700eed61c2e68e37bb0607d4236ada34c602ed86e0160b0146a9e1b339",
	"595fe35df7e10174845e5ce6f334f7be092a82b56bdaccbbf19bb5b11d33d171",

	"78cb0a784a503a4386928f22754e772223587faf2cc8950dd6e386680a86ba65",
	"a4b5c97f1972b54c13d7b050666d6d198ad63e21a9c06e965025ca7e44e9f000",
	"5764e56b2d8896fedf97c30fc2a9d5ad5c2900a0389dd6234c55f166a036ad18",
	"500269aaa0bfe99b6f16a1cd2dd1b21668c5e085b1782de1307449732009cd1b",
	"b936238a475ae0851b89b4c63e094791ac6a03846afad8c7300c77a456bd6e58",
	"6a7690687763248b89a41ac65607a3a88c539a7cfdd3079102b33b2ea607506b",

	"c671c46c3f35165106b26b49873a17a202004883c5a56962d5a13fc4dbcaec44",
	"a04f9846c83204d54f26bab44530c83b4aaa5dfa7772084fcc30a1ab1a9bff1c",
	"552db8f5fad91c4b2d81a021e0e5ef3683c279d8ffa40097417e2fb483c3e23c",
	"aeb8c01dd8a194abb8a6a34c482627e8acbb3b64c8b5470edb9ebddfadc4c61a",
	"8e57d5813af955b94447161bb9965da513f4595e6dc961767cb123657aa1bd50",
	"5b5a5a734813bd45fc5c08015eaca263fa715584b15409730a64b4c20123af7a",
}

// Point multiplication tests. Each group of 3 strings is: P, n, n*P
// Scalar n is not reduced modulo r.
var kat_JQ255S_POINT_MUL = []string {
	"f0f9c911813c90d671dd55050b359c59daa0d39aef85a70196e99c1c1d2bf749",
	"d3a291da69fc56555b79b5f3b6e952650fa69b1168f6d3ea537f78d1043862b5",
	"52c4ef6720b3b137f4a13de2fa9bbed0f88a76b1370ef41860e9cbb50ec1c719",

	"c89aa71a7bbef5a08dd239d6859a59471df19df68a748e98408cdf7b2dc4996c",
	"c7ee606ea14c57307b04ac1d016ad82a2fb1f392c69d15220d701a3d5ed152b4",
	"b961c9ef863fa1fdf8e22dcecb287cb75c8b6bb47cf6ef28f34d12d6daa5b244",

	"71a6c228e378ea8eca58d804aede2eb76f12f4f2f3615963b1e6ff710fc86e32",
	"81fdf930239624fba359fd051ce7f0ae332642331294b911173e83f62299c4fa",
	"053f231139d7f455a3090b12d31f0891929ee11a4e4a9b97b449955708276c46",

	"5f683af7e7132a8c18c7ac224a5018319c646c4482475570d4450dd8ef0e0a7a",
	"adba39cb44b7a97eb45ab0f2b913189d60968d878ff3eb023eefb39b67640bcf",
	"e71c1d259581dd2979f86e054c62f0b5af2b4789fb60685214606c124fd9b90a",

	"5d10110ee82efbc1a71e2368d4adc4ebb019500b26cb643301397f6e6dd7c41c",
	"5a1dd472a60feb189f0390fc83f771847248b7d198a9864e3402cbbf743f7b04",
	"713806b513462b8abbbae0b5adee4f39895b6886e26b6214125502a233f92e17",

	"e117f422cbccf2f7a41a423e8877ee2e6e3e731403e2f39bce089e0b6539bd4e",
	"bdc00054d2123585f4350212829c6677ea5f4a75233621c51137ce0288655c55",
	"502d391351343ab9feb5e6bca688426a06712d91cf1eed186e8beb9de306a41a",

	"fe7a1f889fb102e85cf3ee00bb85d632cdfb6f7c74d112c5e9eb803df80cfc29",
	"ff8f18a9b5fd7ee93676b3a0328aa2a7920db119e6a037247a5e44b0bd52154e",
	"fbab02af2715b10115554136d6d189a8182af4c02585d3c1e40e0bcc214ab700",

	"dcc05a8fa09832a6b63f6f02ebf7d682952b53da2f199ca7ebb8717da81c9f43",
	"de9e8923f07ac1aa7582b4a4c4b4988621f7a5f0aec01f7cb72cf3b089d37d2a",
	"a80183ba294a31c8a42cb6bc48829ebeb53c3279da3990d3597e4fb225472238",

	"a4ed79ad4684a02015a6e68b40b794a076f2b75ab3f1b7149423fc052988c138",
	"977e9951198cad5198c94dba80835fc38d8936532231360b374dcca50f6478c1",
	"f89604a1ee8c0fc6a57efc7d2a6aafb7952defc745ae79316ec964c49a155b16",

	"e02d77b0cd833ff914493d433e66de18368dfc96230541447c67becd22f34a3f",
	"f6f0ad692854a7438d09a0a6c78abff0eaabcff20e02674e9dc2c404960b84c6",
	"3f94fa51d6b6765e97be13a65386a7a7b64321c4e9ced88722f7ac1947ddb902",

	"ef16ce28b6681002c0086d1c1e99b4a7cf2d5bd297901268c8b73491f6607e1c",
	"68c9a4801ea204ae47f7e1dbc8ed47343583c7c9d4c5f1d2c8a9ae7f6e140fae",
	"0b094bd542729485eb7661574f38e2ef758dd1dd0efb9e0b75c41e89bebe5a00",

	"78143f228afd27f6b4b179b32e9eea3ff600ba636515462b83a6d0a23ce7794e",
	"333ff5a622f4565e65f5c070867f1fcae8c75f3de7ad585386dba79465e5eef3",
	"1e1f4ab65f93f824bc0da349a63a467494e2050615aef1d49bafe9402858e864",

	"fd4b6116d3f096fa2922c9e350322558c371ed67c1d1ddccf82918bc9f0db34e",
	"4a461514ff2738927d0138f53123bb9111fa46520c4c83347ad672e72d94be0e",
	"5f0d63130bc98a59dbb7d537ff42c923866f51b529ade1978811ed1f91c8cd60",

	"5367b651b5bc251a732b1314cfff3349187db32563e08516f663950485263d6b",
	"b6c0c6f08e666c9d74a5592f1d47cdf3757549139cbdd4d3ff9c2d2651d68b34",
	"5d73e39320332aaddd578b61c2724d6edd09e4e8e41ca6bfe23dbc796a79f07e",

	"e19574765ff50b7a051b58492b1b1d2ee79b77e3ae9adc4a7be3c1c2c2cbd611",
	"b18f05606bc26ec86adaa37467ec53a87219f2013c45369e8104efaf6766cd13",
	"4eaf39f207d63efbf6346dd326d58c8caba08297037674dc5883d04bbe94d969",

	"2a0b900fdd230363e10ec04b5a744ea51e3535ab378beeb240730b90ea8d9717",
	"88f8d56c3bbbdea6ff663aa9c8a3bd42208c667b7ffe5e46cef9deeb7cd25ca0",
	"dd81f61622406f03b597d8b2a679c63502721bb28fd4e9439f034ee45c171148",

	"2fe4113ca1c62be99012a5f6d55014db36e692974834226952385fceeb37e814",
	"12a60632c5ac702f8ab08ca6f5a64cc016a8b7dfdebb4c333329f7f8674c768e",
	"05cecc583d2812f347939f2edaba68dae21501b5e454c3ea1b104f41f46cd41d",

	"0df7af74b61091a0857017b72f480196acbe1bf2f2d404463989b9ae9a6b7874",
	"36e01111ffc63d5c5605022baba14b95687f1a6d0b09ea793148ab5a58d9e344",
	"153ab43c58ffc4bb25347f65bf422bb617ef67f661e9cc49253969ff3bb10d75",

	"f1a2cbfe94aa32dd9124d88ccc6077441ab119770618642a9414389014ae5a52",
	"678753e1d1de60d2f3622d949e686a7560809737bd806c2fce29bd16c6c11117",
	"9c36fc160d3cfb33508f16ae0e3c746c16640ee9686193d3cfa1581cd55f0243",

	"530471728dc9c1970c472ec35cc542eed4cf87dfb6b422f637ef993cccfe7165",
	"6051a19b837b3717f9bc1e04bf730b87ee0b5ddaaa4d85cb20b4f6f5a15b0665",
	"5f13501b2224e084046778ce9f907cac3393e204584546a1520fdb6cd4e11e6f",
}

var kat_JQ255S_MC_POINT_MUL = []string {
	"a0df5043c5cf6695dd10e3492495821b68457cb2645979ca2fb3c936544d2a18",

	"ba6f94d4e9120982a4a821a3335223d33d89ba373a3a91d4c5cbab3cc28fb252",
	"db19a390cbb3266ca90368abdd318dd056de85f5168766fffe1e54fc35dc9a2a",
	"6c597c2599ab716abded11efb0db0a3683bf7b7ea14ac9dca2a7c883f0340305",
	"f1df069647f98205fe586cd24db4a10bc3965cbed45ae44768072fc530deaf68",
	"b1244425b1281e2e653a8b7edad2b346b616590b01f492f959658ceb50425e0f",
	"4b6085996fe88129101015c3dbfe0614cb53f26e7e5cc35b37635acdabf9ce22",
	"9cc8c0d0fdef61e75f5048f6c2996975a67afe5ba65b93a00683e5cb27ce2a52",
	"8915b6e1727b5e8f566d32637ee1f8c64289bfd0b252a6fd00c2162a4bf36125",
	"0b5596fbc23ae6954c077e3e0c910de0bda71fec86b14493166f0b784e085645",
	"8c5464a02b64b6e467f5f26cf39bc817e9accbe8fe41d8d58bbc5ae1ad8a4f1a",
}
