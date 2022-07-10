# 定义素数域运算 
# 重载了运算符；
# Fp 类型可以与 Fp类型 和 int 类型进行运算；
# 输出结果均为 Fp 类型。

def FiniteField(p):
    class Fp:
        def __init__(self, val: int):
            if(isinstance(val, int)):
                self.val = val
            else:
                self.val = val.val
        
        # self + other  
        def __add__(self, other):    
            if (isinstance(other, int)):
                return Fp((self.val + other) % Fp.p)
            else:
                return Fp((self.val + other.val) % Fp.p)
        
        def __radd__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val + other) % Fp.p)
            else:
                return Fp((self.val + other.val) % Fp.p)
        
        # self - other 
        def __sub__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val - other) % Fp.p)
            else:
                return Fp((self.val - other.val) % Fp.p)
        # other - self
        def __rsub__(self, other):
            if (isinstance(other, int)):
                return Fp((other - self.val) % Fp.p)
            else:
                return Fp((other.val - self.val) % Fp.p)
        
        # * 
        def __mul__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val * other) % Fp.p)
            else:
                return Fp((self.val * other.val) % Fp.p)
        def __rmul__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val * other) % Fp.p)
            else:
                return Fp((self.val * other.val) % Fp.p)
        
        # /
        def __truediv__(self, other):
            if (isinstance(other, int)):
                return self * Fp(other).inv()
            else:
                return self * other.inv()
        def __rtruediv__(self, other):
            if (isinstance(other, int)):
                return self.inv() * Fp(other)
            else:
                return self.inv() * other

        # self ** other       
        def __pow__(self, other):
            if (isinstance(other, int)):
                return Fp(pow(self.val, other, Fp.p))
            else:
                return Fp(pow(self.val, other.val, Fp.p))
        
        # other ** self       
        def __rpow__(self, other):
            if (isinstance(other, int)):
                return Fp(pow(other, self.val, Fp.p))
            else:
                return Fp(pow(other.val, self.val, Fp.p))
        
        # ==
        def __eq__(self, other):
            if (isinstance(other, int)):
                return self.val == other
            else:
                return self.val == other.val
        
        # !=
        def __ne__(self, other):
            if (isinstance(other, int)):
                return not (self.val == other)
            else:
                return not (self.val == other.val)

        # <
        def __lt__(self, other):
            if (isinstance(other, int)):
                return (self.val < other)
            else:
                return (self.val < other.val)
        
        # <=
        def __le__(self, other):
            if (isinstance(other, int)):
                return (self.val <= other)
            else:
                return (self.val <= other.val)
        
        # >
        def __gt__(self, other):
            if (isinstance(other, int)):
                return (self.val > other)
            else:
                return (self.val > other.val)

        # >=
        def __ge__(self, other):
            if (isinstance(other, int)):
                return (self.val >= other)
            else:
                return (self.val >= other.val)
        
        # &
        def __and__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val & other) % Fp.p)
            else:
                return Fp((self.val & other.val) % Fp.p)
        def __rand__(self, other):
            return self.__and__(self, other)
        
        # |
        def __or__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val | other) % Fp.p)
            else:
                return Fp((self.val | other.val) % Fp.p)
        def __ror__(self, other):
            return self.__or__(self, other)
        
        # ^
        def __xor__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val ^ other) % Fp.p)
            else:
                return Fp((self.val ^ other.val) % Fp.p)
        def __rxor__(self, other):
            return self.__xor__(self, other)

        # self >> other
        def __rshift__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val >> other) % Fp.p)
            else:
                return Fp((self.val >> other.val) % Fp.p)
        
        # other >> self
        def __rrshift__(self, other):
            if (isinstance(other, int)):
                return Fp((other >> self.val) % Fp.p)
            else:
                return Fp((other.val >> self.val) % Fp.p)
        
        # self << other
        def __lshift__(self, other):
            if (isinstance(other, int)):
                return Fp((self.val << other) % Fp.p)
            else:
                return Fp((self.val << other.val) % Fp.p)
        
        # other << self
        def __rlshift__(self, other):
            if (isinstance(other, int)):
                return Fp((other << self.val) % Fp.p)
            else:
                return Fp((other.val << self.val) % Fp.p)
        
        # print
        def __repr__(self):
            return hex(self.val)
        
        # int
        def __int__(self):
            return int(self.val)
        
        def inv(self):
        # Field inverse (inverse of 0 is 0).
            return Fp(pow(self.val, Fp.p - 2, Fp.p))
        
        def iszero(self):
            return self.val == 0
        
        def to_bytes(self, length, byteorder):
            return self.val.to_bytes(length, byteorder)
        
        def from_bytes(self, bytes, byteorder):
            self.val = self.val.from_bytes(bytes, byteorder) % Fp.p
    Fp.p = p
    return Fp

if __name__ == "__main__":

    import random

    p = 2**255 - 19
    Fp = FiniteField(p)
    x = Fp(random.randint(0, p-1))
    y = Fp(random.randint(0, p-1))

    t = x.val

    assert (x+y).val == (x+y.val).val
    assert (x*y).val == (x*y.val).val
    assert (x**y).val == (x**y.val).val
    z = x / y
    x >> y

    x << 10

    z1 = y / x
    z2 = y / t

    z3 = x / y
    z4 = t / y

    assert z1 == z2
    assert z3 == z4


