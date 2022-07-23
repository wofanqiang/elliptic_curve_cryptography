import math


print(f"\n------------ 2^(2n) mod p -------------")
def mmm_b(x, y, q, e):
    w = 1
    r = 0
    for i in range(e):
        yi = y & 1
        x0 = x & 1
       # x0 = x
        r0 = r & 1
        u = ((r0+yi*x0)*w)%2
        r = (r + yi*x + u*q)>>1
        y = y >> 1
    if(r>= 2**e): r = r - q
    return r


# little trick
n = 2048
p = 2**n - 17

r = (2**(n-1) + 2**(n-1)) % p

rr = r
x = 16

# rr = r * 2^x
for i in range(x):
    rr = (rr + rr) % p 

# rr = r * 2^(2k)
y = int(math.log2(n) - math.log2(x))
for i in range(y):
    rr = mmm_b(rr, rr, p, n)

print(f"act_rrp = {hex(rr % p)}")
print(f"req_rrp = {hex(2**(2*n) % p)}")
assert rr  == 2**(2*n) % p


#Fast inverse trick -p^(-1) mod 2 ^ k

print(f"\n------------ -p^(-1) mod 2 ^ k -------------")
def get_mont_param(p, bit_len):
    import math
    x = 1
    y = p % (2**bit_len)
    n = int(math.log2(bit_len))

    for i in range(n):
        x = x * (2-y*x) % (2**bit_len)
    
    return (2**bit_len) - x

#p = 0xFFFFFFFE_FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFF_00000000_FFFFFFFF_FFFFFFFF

p = 2**255 - 19

x = get_mont_param(p, 64)

print(f"mont_param = {hex(x)}")
