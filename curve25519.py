from fp import *

class curve25519:
    def __init__(self, u):

        self.p = 2**255 - 19
        self.Fp = FiniteField(self.p)
        self.bits = 255
        self.a24 = self.Fp(121665)
        
        
        if (isinstance(u, int)):
            self.u = u
            self.bytes = self.u.to_bytes(32, 'big')
        elif (isinstance(u, bytes) and len(u)<=32):
            self.bytes = u
            self.u = int.from_bytes(self.bytes, 'big')
        else:
            return False
        
        self.real = self.Fp(self.decodeUCoordinate(self.bytes))
        
    def __repr__(self):
        return hex(self.u)

    def __int__(self):
        return int(self.u)

    def __mul__(self, other):
        if (isinstance(other, int)):
            k = self.decodeScalar(other.to_bytes(32, 'big'))
        elif (isinstance(other, bytes) and len(other)<=32):
            k = self.decodeScalar(other)
        else:
            return False
        res = self.mul(k)
        res = self.encodeUCoordinate(res)
        return curve25519(res)
    
    def __rmul__(self, other):
        return self.__mul__(other)
        
    def decodeLittleEndian(self, b):
        return sum([ b[i] << 8*i for i in range((self.bits+7)//8) ])
    
    def decodeUCoordinate(self, u):
        u_list = [b for b in u]
        # Ignore any unused bits.
        if self.bits % 8:
            u_list[-1] &= (1 << (self.bits % 8)) - 1
        return self.decodeLittleEndian(u_list)
    
    def encodeUCoordinate(self, u):
        #return bytearray([ (u.val >> 8*i) & 0xff for i in range((self.bits+7)//8) ])
        return u.val.to_bytes(32, 'little')

    # k 为 以byte为单位划分的 2进制数据 
    def decodeScalar(self, k):
        k_list = [b for b in k]
        k_list[0] &= 248
        k_list[31] &= 127
        k_list[31] |= 64
        return self.decodeLittleEndian(k_list)
    
    def cswap(self, swap, x_2, x_3):
        "Conditional swap in constant time."
        dummy = swap * (x_2 - x_3)
        x_2 = x_2 - dummy
        x_3 = x_3 + dummy
        return x_2, x_3
    
    def mul(self, k):
        Fp = FiniteField(self.p)
        x_1 = self.real
        x_2 = self.Fp(1)
        z_2 = self.Fp(0)
        x_3 = self.real
        z_3 = self.Fp(1)
        swap = 0

        for t in range(self.bits-1, -1, -1):
            k_t = (k >> t) & 1
            swap ^= k_t
            (x_2, x_3) = self.cswap(swap, x_2, x_3)
            (z_2, z_3) = self.cswap(swap, z_2, z_3)
            swap = k_t

            A = x_2 + z_2
            AA = A**2
            B = x_2 - z_2
            BB = B**2
            E = AA - BB
            C = x_3 + z_3
            D = x_3 - z_3
            DA = D * A
            CB = C * B
            x_3 = (DA + CB)**2
            z_3 = x_1 * (DA - CB)**2
            x_2 = AA * BB
            z_2 = E * (AA + self.a24 * E)

        x_2, x_3 = self.cswap(swap, x_2, x_3)
        z_2, z_3 = self.cswap(swap, z_2, z_3)
        res = x_2 / z_2
        return res

if __name__ == "__main__":
    u = curve25519(0xe6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c)
    k = 0xa546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4.to_bytes(32, 'big')
    u_new = k*u
    print(u_new.real)
    print(u.real)