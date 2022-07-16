from inspect import signature
from fp import *
import hashlib


class ed448_point:
    def __init__(self, s):
        self.p = 2**448 - 2**224 - 1
        self.Fp = FiniteField(self.p)
        self.d = self.Fp(-39081)
        if(isinstance(s, int)):
            self.s_bytes = s.to_bytes(57, 'big')
            self.s_int = s
        elif(isinstance(s, bytes)):
            self.s_bytes = s
            self.s_int = int.from_bytes(s, 'big')
        else:
            raise TypeError("Only Int or Bytes.")
        
        (self.X, self.Y, self.Z) = self.point_decompress(self.s_bytes)
        self.decode_fail = (self.X is None)
        self.x = self.X
        self.y = self.Y
        self.compressed_val = s
    
    def __eq__(self, other):
        if(isinstance(other, ed448_point)):
            raise TypeError("Not a point.")
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (self.X * other.Z - other.X * self.Z  != 0):
            return False
        if (self.Y * other.Z - other.Y * self.Z != 0) :
            return False
        return True
    
    def __neq__(self, other):
        if(isinstance(other, ed448_point)):
            return True
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (self.X * other.Z - other.X * self.Z  == 0 and self.Y * other.Z - other.Y * self.Z == 0):
            return False
        
        return True
    
    def point_decompress(self, s):
        if len(s) != 57:
            raise Exception("Invalid input length for decompression")
        y = int.from_bytes(s, "little")
        sign = y >> 455
        y = self.Fp(y & (1 << 448) - 1)

        if(y >= self.p):
            return (None, None, None)
        
        x = self.recover_x(y, sign)
        
        if x is None:
            return (None, None, None)
        else:
            return (x, y, 1)

    def point_compress(self):
        x = self.x.val
        y = self.y.val
        return int.to_bytes(y | ((x & 1) << 455), 57, "little")
    
    def recover_x(self, y, sign):
        if y >= self.p:
            return None
        
        u = y*y -1
        v = self.d*y*y-1
        x = (u**3 * v) * ((u**5 * v**3) ** ((self.p-3) // 4))
        
        if (v*x*x - u) == 0:
            x = x
        else:
            return None
        
        if(x == 0 and sign == 1):
            return None
        elif (x & 1) != sign:
            x = self.p - x

        return x

class ed448_point_xy(ed448_point):
    def __init__(self, x, y):
        self.p = 2**448 - 2**224 - 1
        self.Fp = FiniteField(self.p)
        self.x = self.Fp(x)
        self.y = self.Fp(y)
        self.d = self.Fp(-39081)
        # Square root of -1
        self.X, self.Y, self.Z = self.x, self.y, self.Fp(1)
        self.compressed_val = self.point_compress()

class ed448_point_xyz(ed448_point):
    def __init__(self, X, Y, Z):

        self.p = 2**448 - 2**224 - 1
        self.Fp = FiniteField(self.p)
        self.X, self.Y, self.Z = self.Fp(X), self.Fp(Y), self.Fp(Z)
        self.x = self.X / self.Z
        self.y = self.Y / self.Z
        self.d = self.Fp(-39081)
        # Square root of -1
        self.modp_sqrt_m1 = self.Fp(2) ** ((self.p-1) // 4)
        self.compressed_val = self.point_compress()

class ed448:
    def __init__(self):
        # Base field Z_p
        self.p = 2**448 - 2**224 - 1
        self.Fp = FiniteField(self.p)
        # Curve constant
        self.d = self.Fp(-39081)
        self.n = 447
        self.a = self.Fp(1)
        self.b = 456
        self.c = 2
        self.Bx = 224580040295924300187604334099896036246789641632564134246125461686950415467406032909029192869357953282578032075146446173674602635247710
        self.By = 298819210078481492676017930443930673437544040154080242095928241372331506189835876003536878655418784733982303233503462500531545062832660
        self.G = ed448_point_xy(self.Bx, self.By)
        # Group order
        self.q = 2**446 - 13818066809895115352007386748515426880336692474882178609894547503885
        

    def shake256(self, s, dlen):
        ph = hashlib.shake_256()
        ph.update(s)
        return ph.digest(dlen)
    
    def shake256_modq(self, s, dlen):
        return int.from_bytes(self.shake256(s, dlen), "little") % self.q

    # Points are represented as tuples (X, Y, Z, T) of extended
    # coordinates, with x = X/Z, y = Y/Z, x*y = T/Z
    def point_add(self, P: ed448_point_xyz, Q: ed448_point_xyz):
        A = P.Z * Q.Z
        B = A**2
        C = P.X * Q.X
        D = P.Y * Q.Y
        E = self.d * C * D
        F = B - E
        G = B + E
        H = (P.X + P.Y) * (Q.X + Q.Y)
        X3 = A * F * (H - C - D)
        Y3 = A * G * (D - C)
        Z3 = F * G
        return ed448_point_xyz(X3, Y3, Z3)
    def point_double(self, P: ed448_point_xyz):
        B = (P.X + P.Y) ** 2
        C = P.X ** 2
        D = P.Y ** 2
        E = C + D
        H = P.Z ** 2
        J = E - 2 * H
        X3 = (B - E) * J
        Y3 = E * (C - D)
        Z3 = E * J
        return ed448_point_xyz(X3, Y3, Z3)
    
    def point_double_xyz(self, P: ed448_point_xyz):
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
        return ed448_point_xyz(X3, Y3, Z3)  

    def point_add_xyz(self, P: ed448_point_xyz, Q: ed448_point_xyz):
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
        return ed448_point_xyz(X3, Y3, Z3) 

    # Computes Q = s * Q
    def point_mul(self, s, P: ed448_point):

        Q = ed448_point_xyz(0, 1, 1)  # Neutral element

        while s > 0:
            if s & 1:
                Q = self.point_add(Q, P)
            P = self.point_double(P)
            s >>= 1
        return Q
    
    def point_equal(sel, P: ed448_point, Q: ed448_point):
        # x1 / z1 == x2 / z2  <==>  x1 * z2 == x2 * z1
        if (P.X * Q.Z - Q.X * P.Z  != 0):
            return False
        if (P.Y * Q.Z - Q.Y * P.Z != 0) :
            return False
        return True
    
    def secret_expand(self, secret):
        if len(secret) != 57:
            raise Exception("Bad size of private key")
        h = self.shake256(secret, 114)
        a = int.from_bytes(h[:57], "little")
        a &=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc
        #TODO: the highest bit of the second to last octet is set. 进行该运算会出错
        #a |= 0x8080808080808080808080808080808080808080808080808080808080808080808080808080808080808080808080808080808080808000
        return (a, h[57:])
    
    def generate_public_key(self, secret):
        (a, dummy) = self.secret_expand(secret)
        pub_key = self.point_mul(a, self.G)
        return pub_key.compressed_val
    
    def dom4(self, x: int, y: bytes):
        #"SigEd448" is in ASCII (8 octets)
        const_val = 0x5369674564343438
        const_val = const_val.to_bytes(8, 'big')
        x_byte = x.to_bytes(1, 'big')
        y_len = len(y).to_bytes(1, 'big')
        #return "SigEd448"|| octet(x) || octet(OLEN(y)) || y
        return const_val + x_byte + y_len + y

    def sign(self, secret: bytes, msg: bytes, phflag=0, context = b''):
        #Compute public key
        a, prefix = self.secret_expand(secret)
        A = self.point_mul(a, self.G).compressed_val
        #compute sign
        if(phflag == 0):
            r = self.shake256_modq(self.dom4(phflag, context) + prefix + msg, 114)
        elif(phflag == 1):
            r = self.shake256_modq(self.dom4(phflag, context) + prefix + self.shake256(msg, 64), 114)
        else:
            return False

        R = self.point_mul(r, self.G)
        Rs = R.compressed_val
        if(phflag == 0):
            k = self.shake256_modq(self.dom4(phflag, context) + Rs + A + msg, 114)
        elif(phflag == 1):
            k = self.shake256_modq(self.dom4(phflag, context) + Rs + A + self.shake256(msg, 64), 114)
        else:
            return False
        
        s = (r + k * a) % self.q
        return Rs + int.to_bytes(s, 57, "little")
    
    
    def verify(self, public, msg, signature, phflag=0, context = b''):
        if len(public) != 57:
            raise Exception("Bad public key length")
        if len(signature) != 114:
            raise Exception("Bad signature length")
        #Decode the public key as point A (point_decompress())
        A = ed448_point(public)
        #If the decoding fail the signature is invalid
        if A.decode_fail:
            return False

        #Decode the Rs as point R (point_decompress())
        Rs = signature[:57]
        R = ed448_point(Rs)
        #If the decoding fail the signature is invalid
        if R.decode_fail:
            return False

        #Check the signature s 
        s = int.from_bytes(signature[57:], "little")
        if s >= self.q: 
            return False
        
        if(phflag == 0):
            k = self.shake256_modq(self.dom4(phflag, context) + Rs + public + msg, 114)
        elif(phflag == 1):
            k = self.shake256_modq(self.dom4(phflag, context) + Rs + public + self.shake256(msg, 64), 114)
        else:
            return False

        sB = self.point_mul(s, self.G)
        kA = self.point_mul(k, A)
        return self.point_equal(sB, self.point_add(R, kA))


if __name__ == "__main__":
    ed = ed448()


    print("\nRFC8032 ED448 TEST 1")
    sk = 0x6c82a562cb808d10d632be89c8513ebf6c929f34ddfa8c9f63c9960ef6e348a3528c8a3fcc2f044e39a3fc5b94492f8f032e7549a20098f95b
    sk = sk.to_bytes(57, 'big')

    req_pk = 0x5fd7449b59b461fd2ce787ec616ad46a1da1342485a70e1f8a0ea75d80e96778edf124769b46c7061bd6783df1e50f6cd1fa1abeafe8256180
    req_pk = req_pk.to_bytes(57, 'big')
    req_sign = 0x533a37f6bbe457251f023c0d88f976ae2dfb504a843e34d2074fd823d41a591f2b233f034f628281f2fd7a22ddd47d7828c59bd0a21bfd3980ff0d2028d4b18a9df63e006c5d1c2d345b925d8dc00b4104852db99ac5c7cdda8530a113a0f4dbb61149f05a7363268c71d95808ff2e652600
    req_sign = req_sign.to_bytes(114, 'big')

    pk = ed.generate_public_key(sk)

    assert req_pk == pk
    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    
    msg = b''
    print(f"msg = {msg.hex()}")

    phflag=0
    act_sign = ed.sign(sk, msg, phflag)
    assert req_sign == act_sign
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag)
    print(f"Verify: {verify_status}")  



    
    print("\nRFC8032 ED448 TEST 2")
    sk = 0xc4eab05d357007c632f3dbb48489924d552b08fe0c353a0d4a1f00acda2c463afbea67c5e8d2877c5e3bc397a659949ef8021e954e0a12274e
    sk = sk.to_bytes(57, 'big')

    req_pk = 0x43ba28f430cdff456ae531545f7ecd0ac834a55d9358c0372bfa0c6c6798c0866aea01eb00742802b8438ea4cb82169c235160627b4c3a9480
    req_pk = req_pk.to_bytes(57, 'big')
    req_sign = 0xd4f8f6131770dd46f40867d6fd5d5055de43541f8c5e35abbcd001b32a89f7d2151f7647f11d8ca2ae279fb842d607217fce6e042f6815ea000c85741de5c8da1144a6a1aba7f96de42505d7a7298524fda538fccbbb754f578c1cad10d54d0d5428407e85dcbc98a49155c13764e66c3c00
    req_sign = req_sign.to_bytes(114, 'big')

    pk = ed.generate_public_key(sk)

    assert req_pk == pk
    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    
    msg = 0x03
    msg = msg.to_bytes(1, 'big')
    print(f"msg = {msg.hex()}")

    context = 0x666f6f
    context = context.to_bytes(3, 'big')

    phflag=0
    act_sign = ed.sign(sk, msg, phflag, context)
    assert req_sign == act_sign
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag, context)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag, context)
    print(f"Verify: {verify_status}")  


    print("\nRFC8032 ED448hp TEST 1")
    sk = 0x833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42ef7822e0d5104127dc05d6dbefde69e3ab2cec7c867c6e2c49
    sk = sk.to_bytes(57, 'big')

    req_pk = 0x259b71c19f83ef77a7abd26524cbdb3161b590a48f7d17de3ee0ba9c52beb743c09428a131d6b1b57303d90d8132c276d5ed3d5d01c0f53880
    req_pk = req_pk.to_bytes(57, 'big')
    req_sign = 0x822f6901f7480f3d5f562c592994d9693602875614483256505600bbc281ae381f54d6bce2ea911574932f52a4e6cadd78769375ec3ffd1b801a0d9b3f4030cd433964b6457ea39476511214f97469b57dd32dbc560a9a94d00bff07620464a3ad203df7dc7ce360c3cd3696d9d9fab90f00
    req_sign = req_sign.to_bytes(114, 'big')

    pk = ed.generate_public_key(sk)

    assert req_pk == pk
    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    # 'abc'
    msg = 0x616263
    msg = msg.to_bytes(3, 'big')
    print(f"msg = {msg.hex()}")

    context = 0x666f6f
    context = context.to_bytes(3, 'big')

    phflag=1
    act_sign = ed.sign(sk, msg, phflag)
    assert req_sign == act_sign
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag)
    print(f"Verify: {verify_status}")  



    print("\nRFC8032 ED448hp TEST 2")
    sk = 0x833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42ef7822e0d5104127dc05d6dbefde69e3ab2cec7c867c6e2c49
    sk = sk.to_bytes(57, 'big')

    req_pk = 0x259b71c19f83ef77a7abd26524cbdb3161b590a48f7d17de3ee0ba9c52beb743c09428a131d6b1b57303d90d8132c276d5ed3d5d01c0f53880
    req_pk = req_pk.to_bytes(57, 'big')
    req_sign = 0xc32299d46ec8ff02b54540982814dce9a05812f81962b649d528095916a2aa481065b1580423ef927ecf0af5888f90da0f6a9a85ad5dc3f280d91224ba9911a3653d00e484e2ce232521481c8658df304bb7745a73514cdb9bf3e15784ab71284f8d0704a608c54a6b62d97beb511d132100
    req_sign = req_sign.to_bytes(114, 'big')

    pk = ed.generate_public_key(sk)

    assert req_pk == pk
    print(f"sk = {sk.hex()}")
    print(f"pk = {pk.hex()}")

    # 'abc'
    msg = 0x616263
    msg = msg.to_bytes(3, 'big')
    print(f"msg = {msg.hex()}")

    context = 0x666f6f
    context = context.to_bytes(3, 'big')

    phflag=1
    act_sign = ed.sign(sk, msg, phflag, context)
    assert req_sign == act_sign
    print(f"sign = {act_sign.hex()}")

    verify_status = ed.verify(pk, msg, act_sign, phflag, context)
    print(f"Verify: {verify_status}")

    verify_status = ed.verify(pk, msg+b'1', act_sign, phflag, context)
    print(f"Verify: {verify_status}")  