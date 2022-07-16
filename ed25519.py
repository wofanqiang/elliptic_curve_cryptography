from inspect import signature
from fp import *
import hashlib


class ed25519_point:
    def __init__(self, s):
        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        self.d = self.Fp(-121665) / self.Fp(121666)
        self.modp_sqrt_m1 = self.Fp(2) ** ((self.p-1) // 4)
        if(isinstance(s, int)):
            self.s_bytes = s.to_bytes(32, 'big')
            self.s_int = s
        elif(isinstance(s, bytes)):
            self.s_bytes = s
            self.s_int = int.from_bytes(s, 'big')
        else:
            raise TypeError("Only Int or Bytes.")
        
        (self.X, self.Y, self.Z, self.T) = self.point_decompress(self.s_bytes)
        self.decode_fail = (self.X is None)
        self.x = self.X
        self.y = self.Y
        self.compressed_val = s
    
    def __eq__(self, other):
        if(isinstance(other, ed25519_point)):
            raise TypeError("Not a point.")
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (self.X * other.Z - other.X * self.Z  != 0):
            return False
        if (self.Y * other.Z - other.Y * self.Z != 0) :
            return False
        return True
    
    def __neq__(self, other):
        if(isinstance(other, ed25519_point)):
            return True
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (self.X * other.Z - other.X * self.Z  == 0 and self.Y * other.Z - other.Y * self.Z == 0):
            return False
        
        return True
    
    def point_decompress(self, s):
        if len(s) != 32:
            raise Exception("Invalid input length for decompression")
        y = int.from_bytes(s, "little")
        sign = y >> 255
        y = self.Fp(y & (1 << 255) - 1)
        
        x = self.recover_x(y, sign)
        
        if x is None:
            return (None, None, None, None)
        else:
            return (x, y, 1, x*y)

    def point_compress(self):
        x = self.x.val
        y = self.y.val
        return int.to_bytes(y | ((x & 1) << 255), 32, "little")
    
    def recover_x(self, y, sign):
        if y >= self.p:
            return None
        
        u = y*y -1
        v = self.d*y*y+1
        x = (u * v**3) * ((u * v**7) ** ((self.p-5) // 8))
        
        if (v*x*x - u) == 0:
            x = x
        elif (v*x*x + u) == 0:
            x = x * self.modp_sqrt_m1
        else:
            return None
        
        if(x == 0 and sign == 1):
            return None
        elif (x & 1) != sign:
            x = self.p - x

        return x

class ed25519_point_xy(ed25519_point):
    def __init__(self, x, y):
        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        self.x = self.Fp(x)
        self.y = self.Fp(y)
        self.d = self.Fp(-121665) / self.Fp(121666)
        # Square root of -1
        self.modp_sqrt_m1 = self.Fp(2) ** ((self.p-1) // 4)
        self.X, self.Y, self.Z, self.T = self.x, self.y, self.Fp(1), self.x*self.y
        self.compressed_val = self.point_compress()

class ed25519_point_xyz(ed25519_point):
    def __init__(self, X, Y, Z):

        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        self.X, self.Y, self.Z = self.Fp(X), self.Fp(Y), self.Fp(Z)
        self.x = self.X / self.Z
        self.y = self.Y / self.Z
        self.d = self.Fp(-121665) / self.Fp(121666)
        # Square root of -1
        self.modp_sqrt_m1 = self.Fp(2) ** ((self.p-1) // 4)
        self.compressed_val = self.point_compress()

class ed25519_point_xyzt(ed25519_point):
    def __init__(self, X, Y, Z, T):

        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        self.X, self.Y, self.Z, self.T = self.Fp(X), self.Fp(Y), self.Fp(Z), self.Fp(T)
        self.x = self.X / self.Z
        self.y = self.Y / self.Z
        self.d = self.Fp(-121665) / self.Fp(121666)
        # Square root of -1
        self.modp_sqrt_m1 = self.Fp(2) ** ((self.p-1) // 4)
        self.compressed_val = self.point_compress()

class ed25519:
    def __init__(self):
        # Base field Z_p
        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        # Curve constant
        self.d = self.Fp(-121665) / self.Fp(121666)
        self.n = 254
        self.a = self.Fp(-1)
        self.b = 256
        self.c = 3
        self.Bx = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        self.By = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        self.G = ed25519_point_xy(self.Bx, self.By)
        # Group order
        self.q = 2**252 + 27742317777372353535851937790883648493
        self.Fq = FiniteField(self.q)

    def sha512(self, s):
        return hashlib.sha512(s).digest()
    
    def sha512_modq(self, s):
        return int.from_bytes(self.sha512(s), "little") % self.q

    # Points are represented as tuples (X, Y, Z, T) of extended
    # coordinates, with x = X/Z, y = Y/Z, x*y = T/Z
    def point_add(self, P: ed25519_point_xyzt, Q: ed25519_point_xyzt):
        A = (P.Y-P.X) * (Q.Y-Q.X)
        B = (P.Y+P.X) * (Q.Y+Q.X)
        C = 2 * P.T * self.d * Q.T
        D = 2 * P.Z * Q.Z
        E = B-A
        F = D-C
        G = D+C
        H = B+A
        X3 = E*F
        Y3 = G*H
        Z3 = F*G
        T3 = E*H
        return ed25519_point_xyzt(X3, Y3, Z3, T3)
    def point_double(self, P: ed25519_point_xyzt):
        A = P.X**2
        B = P.Y**2
        C = 2*P.Z**2
        H = A+B
        E = H-(P.X+P.Y)**2
        G = A-B
        F = C+G
        X3 = E*F
        Y3 = G*H
        T3 = E*H
        Z3 = F*G
        return ed25519_point_xyzt(X3, Y3, Z3, T3)

    def point_double_xyz(self, P: ed25519_point_xyz):
        B = (P.X+P.Y)**2
        C = P.X**2
        D = P.Y**2
        E = self.a*C
        F = E+D
        H = P.Z**2
        J = F-2*H
        X3 = (B-C-D)*J
        Y3 = F*(E-D)
        Z3 = F*J
        return ed25519_point_xyz(X3, Y3, Z3)  

    def point_add_xyz(self, P: ed25519_point_xyz, Q: ed25519_point_xyz):
        A = P.Z*Q.Z
        B = A**2
        C = P.X*Q.X
        D = P.Y*Q.Y
        E = self.d*C*D
        F = B-E
        G = B+E
        X3 = A*F*((P.X+P.Y)*(Q.X+Q.Y)-C-D)
        Y3 = A*G*(D-self.a*C)
        Z3 = F*G
        return ed25519_point_xyz(X3, Y3, Z3)   

    # Computes Q = s * Q
    def point_mul(self, s, P: ed25519_point):

        Q = ed25519_point_xyzt(0, 1, 1, 0)  # Neutral element
        #Q = ed25519_point_xyz(0, 1, 1)  # Neutral element

        while s > 0:
            if s & 1:
                Q = self.point_add(Q, P)
            P = self.point_double(P)
            #    Q = self.point_add_xyz(Q, P)
            #P = self.point_double_xyz(P)
            s >>= 1
        return Q
    
    def point_equal(sel, P: ed25519_point, Q: ed25519_point):
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (P.X * Q.Z - Q.X * P.Z  != 0):
            return False
        if (P.Y * Q.Z - Q.Y * P.Z != 0) :
            return False
        return True
    
    def secret_expand(self, secret):
        if len(secret) != 32:
            raise Exception("Bad size of private key")
        h = self.sha512(secret)
        a = int.from_bytes(h[:32], "little")
        #a &= (1 << 254) - 8
        #a |= (1 << 254)
        a &= 0x3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff8
        a |= 0x4000000000000000000000000000000000000000000000000000000000000000
        return (a, h[32:])
    
    def generate_public_key(self, secret):
        (a, dummy) = self.secret_expand(secret)
        #a = int.from_bytes(a, 'big')
        pub_key = self.point_mul(a, self.G)
        return pub_key.compressed_val
    
    def dom2(self, x: int, y: bytes):
        #"SigEd25519 no Ed25519 collisions" is in ASCII (32 octets)
        const_val = 0x53696745643235353139206e6f204564323535313920636f6c6c6973696f6e73
        const_val = const_val.to_bytes(32, 'big')
        x_byte = x.to_bytes(1, 'big')
        y_len = len(y).to_bytes(1, 'big')
        #return "SigEd25519 no Ed25519 collisions" || octet(x) || octet(OLEN(y)) || y
        return const_val + x_byte + y_len + y

    def sign(self, secret: bytes, msg: bytes, phflag=None, context = b''):
        #Compute public key
        a, prefix = self.secret_expand(secret)
        A = self.point_mul(a, self.G).compressed_val
        #compute sign
        if(phflag == 0):
            r = self.sha512_modq(self.dom2(phflag, context) + prefix + msg)
        elif(phflag == 1):
            r = self.sha512_modq(self.dom2(phflag, context) + prefix + self.sha512(msg))
        else:
            r = self.sha512_modq(prefix + msg)

        R = self.point_mul(r, self.G)
        Rs = R.compressed_val
        if(phflag == 0):
            k = self.sha512_modq(self.dom2(phflag, context) + Rs + A + msg)
        elif(phflag == 1):
            k = self.sha512_modq(self.dom2(phflag, context) + Rs + A + self.sha512(msg))
        else:
            k = self.sha512_modq(Rs + A + msg)
        
        s = (r + k * a) % self.q
        return Rs + int.to_bytes(s, 32, "little")
    
    
    def verify(self, public, msg, signature, phflag=None, context = b''):
        if len(public) != 32:
            raise Exception("Bad public key length")
        if len(signature) != 64:
            raise Exception("Bad signature length")
        #Decode the public key as point A (point_decompress())
        A = ed25519_point(public)
        #If the decoding fail the signature is invalid
        if A.decode_fail:
            return False

        #Decode the Rs as point R (point_decompress())
        Rs = signature[:32]
        R = ed25519_point(Rs)
        #If the decoding fail the signature is invalid
        if R.decode_fail:
            return False

        #Check the signature s 
        s = int.from_bytes(signature[32:], "little")
        if s >= self.q: 
            return False
        
        if(phflag == 0):
            k = self.sha512_modq(self.dom2(phflag, context) + Rs + public + msg)
        elif(phflag == 1):
            k = self.sha512_modq(self.dom2(phflag, context) + Rs + public + self.sha512(msg))
        else:
            k = self.sha512_modq(Rs + public + msg)

        sB = self.point_mul(s, self.G)
        kA = self.point_mul(k, A)
        return self.point_equal(sB, self.point_add_xyz(R, kA))


if __name__ == "__main__":
    ed = ed25519()


    print("\nRFC8032 ED25519 TEST 1")
    sk = 0x9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60
    sk = sk.to_bytes(32, 'big')

    req_pk = 0xd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0xe5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = b''
    print(f"msg = {msg.hex()}")

    act_sign = ed.sign(sk, msg)
    assert req_sign == act_sign

    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign)
    print(f"Verify: {verify_status}")


    print("\nRFC8032 ED25519 TEST 2")
    sk = 0x4ccd089b28ff96da9db6c346ec114e0f5b8a319f35aba624da8cf6ed4fb8a6fb
    sk = sk.to_bytes(32, 'big')

    req_pk = 0x3d4017c3e843895a92b70aa74d1b7ebc9c982ccf2ec4968cc0cd55f12af4660c
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0x92a009a9f0d4cab8720e820b5f642540a2b27b5416503f8fb3762223ebdb69da085ac1e43e15996e458f3613d0f11d8c387b2eaeb4302aeeb00d291612bb0c00
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = 0x72
    msg = msg.to_bytes(1, 'big')
    print(f"msg = {msg.hex()}")

    act_sign = ed.sign(sk, msg)
    assert req_sign == act_sign
    
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign)
    print(f"Verify: {verify_status}")

    print("\nRFC8032 ED25519 TEST 1024")
    sk = 0xf5e5767cf153319517630f226876b86c8160cc583bc013744c6bf255f5cc0ee5
    sk = sk.to_bytes(32, 'big')

    req_pk = 0x278117fc144c72340f67d0f2316e8386ceffbf2b2428c9c51fef7c597f1d426e
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0x0aab4c900501b3e24d7cdf4663326a3a87df5e4843b2cbdb67cbf6e460fec350aa5371b1508f9f4528ecea23c436d94b5e8fcd4f681e30a6ac00a9704a188a03
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = 0x08b8b2b733424243760fe426a4b54908632110a66c2f6591eabd3345e3e4eb98fa6e264bf09efe12ee50f8f54e9f77b1e355f6c50544e23fb1433ddf73be84d879de7c0046dc4996d9e773f4bc9efe5738829adb26c81b37c93a1b270b20329d658675fc6ea534e0810a4432826bf58c941efb65d57a338bbd2e26640f89ffbc1a858efcb8550ee3a5e1998bd177e93a7363c344fe6b199ee5d02e82d522c4feba15452f80288a821a579116ec6dad2b3b310da903401aa62100ab5d1a36553e06203b33890cc9b832f79ef80560ccb9a39ce767967ed628c6ad573cb116dbefefd75499da96bd68a8a97b928a8bbc103b6621fcde2beca1231d206be6cd9ec7aff6f6c94fcd7204ed3455c68c83f4a41da4af2b74ef5c53f1d8ac70bdcb7ed185ce81bd84359d44254d95629e9855a94a7c1958d1f8ada5d0532ed8a5aa3fb2d17ba70eb6248e594e1a2297acbbb39d502f1a8c6eb6f1ce22b3de1a1f40cc24554119a831a9aad6079cad88425de6bde1a9187ebb6092cf67bf2b13fd65f27088d78b7e883c8759d2c4f5c65adb7553878ad575f9fad878e80a0c9ba63bcbcc2732e69485bbc9c90bfbd62481d9089beccf80cfe2df16a2cf65bd92dd597b0707e0917af48bbb75fed413d238f5555a7a569d80c3414a8d0859dc65a46128bab27af87a71314f318c782b23ebfe808b82b0ce26401d2e22f04d83d1255dc51addd3b75a2b1ae0784504df543af8969be3ea7082ff7fc9888c144da2af58429ec96031dbcad3dad9af0dcbaaaf268cb8fcffead94f3c7ca495e056a9b47acdb751fb73e666c6c655ade8297297d07ad1ba5e43f1bca32301651339e22904cc8c42f58c30c04aafdb038dda0847dd988dcda6f3bfd15c4b4c4525004aa06eeff8ca61783aacec57fb3d1f92b0fe2fd1a85f6724517b65e614ad6808d6f6ee34dff7310fdc82aebfd904b01e1dc54b2927094b2db68d6f903b68401adebf5a7e08d78ff4ef5d63653a65040cf9bfd4aca7984a74d37145986780fc0b16ac451649de6188a7dbdf191f64b5fc5e2ab47b57f7f7276cd419c17a3ca8e1b939ae49e488acba6b965610b5480109c8b17b80e1b7b750dfc7598d5d5011fd2dcc5600a32ef5b52a1ecc820e308aa342721aac0943bf6686b64b2579376504ccc493d97e6aed3fb0f9cd71a43dd497f01f17c0e2cb3797aa2a2f256656168e6c496afc5fb93246f6b1116398a346f1a641f3b041e989f7914f90cc2c7fff357876e506b50d334ba77c225bc307ba537152f3f1610e4eafe595f6d9d90d11faa933a15ef1369546868a7f3a45a96768d40fd9d03412c091c6315cf4fde7cb68606937380db2eaaa707b4c4185c32eddcdd306705e4dc1ffc872eeee475a64dfac86aba41c0618983f8741c5ef68d3a101e8a3b8cac60c905c15fc910840b94c00a0b9d0

    msg = msg.to_bytes(1023, 'big')
    print(f"msg = {msg.hex()}")

    act_sign = ed.sign(sk, msg)
    assert req_sign == act_sign
    
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign)
    print(f"Verify: {verify_status}")


    print("\nRFC8032 ED25519 TEST SHA(abc)")
    sk = 0x833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42
    sk = sk.to_bytes(32, 'big')

    req_pk = 0xec172b93ad5e563bf4932c70e1245034c35467ef2efd4d64ebf819683467e2bf
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0xdc2a4459e7369633a52b1bf277839a00201009a3efbf3ecb69bea2186c26b58909351fc9ac90b3ecfdfbc7c66431e0303dca179c138ac17ad9bef1177331a704
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = 0xddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f

    msg = msg.to_bytes(64, 'big')
    print(f"msg = {msg.hex()}")

    act_sign = ed.sign(sk, msg)
    assert req_sign == act_sign
    
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign)
    print(f"Verify: {verify_status}")



    print("\nRFC8032 ED25519 TEST Ed25519ctx")
    sk = 0x0305334e381af78f141cb666f6199f57bc3495335a256a95bd2a55bf546663f6
    sk = sk.to_bytes(32, 'big')

    req_pk = 0xdfc9425e4f968f7f0c29f0259cf5f9aed6851c2bb4ad8bfb860cfee0ab248292
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0x55a4cc2f70a54e04288c5f4cd1e45a7bb520b36292911876cada7323198dd87a8b36950b95130022907a7fb7c4e9b2d5f6cca685a587b4b21f4b888e4e7edb0d
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = 0xf726936d19c800494e3fdaff20b276a8

    msg = msg.to_bytes(16, 'big')
    print(f"msg = {msg.hex()}")

    phflag = 0
    context = 0x666f6f
    context = context.to_bytes(3, 'big')

    act_sign = ed.sign(sk, msg, phflag, context)
    assert req_sign == act_sign
    
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag, context)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag, context)
    print(f"Verify: {verify_status}")


    print("\nRFC8032 ED25519 TEST Ed25519ph")
    sk = 0x833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42

    sk = sk.to_bytes(32, 'big')

    req_pk = 0xec172b93ad5e563bf4932c70e1245034c35467ef2efd4d64ebf819683467e2bf
    req_pk = req_pk.to_bytes(32, 'big')
    req_sign = 0x98a70222f0b8121aa9d30f813d683f809e462b469c7ff87639499bb94e6dae4131f85042463c2a355a2003d062adf5aaa10b8c61e636062aaad11c2a26083406
    req_sign = req_sign.to_bytes(64, 'big')

    pk = ed.generate_public_key(sk)
    assert req_pk == pk

    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    msg = 0x616263
    msg = msg.to_bytes(3, 'big')
    print(f"msg = {msg.hex()}")

    phflag = 1
    context = b''

    act_sign = ed.sign(sk, msg, phflag, context)
    assert req_sign == act_sign
    
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag, context)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag, context)
    print(f"Verify: {verify_status}")