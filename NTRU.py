from poly import Poly
from Zmod import Zmod
from common import *
import time
import random
from tqdm import tqdm
from PIL import Image
import binascii

class NTRU:
    def __init__(self, N, p, q, df, dg, d):
        self.N = N
        self.Rp = p
        self.Rq = q
        self.d = d
        self.df = df
        self.dg = dg
        self.h = None
        self.f = None #private key
        self.g = None #private key
        self.Fp = None
        self.Fq = None
        self.test = None

    def key_gen(self):
        while True:
            try:
                self.f = random_ternary_list(self.df + 1, self.df, self.N - 1)
                self.Fp = Poly(self.Rp, self.f).inv_mod_xN_prime_pow(self.N)
                self.Fq = Poly(self.Rq, self.f).inv_mod_xN_prime_pow(self.N)
                
                # If either inversion fails (returns None), generate a new f and retry
                if self.Fp is not None and self.Fq is not None:
                    break

            except Exception as e:
                continue

        self.g = random_ternary_list(self.dg, self.dg, self.N - 1)
        self.h = Poly(self.Rq, self.g).mul_mod_convolution_ring(self.Fq, self.N)

        return self.h
    
    def encrypt_block(self, f_m: Poly) -> Poly:
        r = random_ternary_list(self.d, self.d, self.N - 1)

        e = (Poly(self.Rq, r) * self.Rp).mul_mod_convolution_ring(self.h, self.N).add_inplace(f_m.change_ring(self.Rq))

        return e
    
    def decrypt_block(self, e: Poly) -> Poly:
        temp = Poly(self.Rq, self.f.copy()).mul_mod_convolution_ring(e, self.N) 
        
        center = temp.center_lift()
        f_m = Poly(self.Rp, center).mul_mod_convolution_ring(self.Fp, self.N)

        return f_m

    # def encrypt(self, m_list: np.ndarray[np.complex128]) -> np.ndarray[np.complex128]:
        
    #     rounding = np.ceil(self.N / 2)
       
    #     pad = lambda x: np.pad(x, pad_width=(0, max(0,int(rounding -  len(x)))), mode='constant', constant_values=0)

    #     polys = complex_array_2_polys(m_list, self.Rp, self.N // 2)
    #     enc_polys = []
    #     for p in polys:
    #         tmp = self.encrypt_block(p)
    #         enc_polys.append(tmp)
    #     self.test = polys
        
    #     enc_polys = polys_2_complex_array(enc_polys,pad)
       
    #     return enc_polys


    # def decrypt(self, enc_list: np.ndarray[np.complex128]) -> np.ndarray[np.complex128]:
        
    #     polys = complex_array_2_polys(enc_list, self.Rq, self.N // 2 + 1)
    #     m_list = [self.decrypt_block(p) for p in polys]
        
    #     return polys_2_complex_array(m_list)
    def encrypt(self, m: list, block_size = 96):
        '''
        m: list of bit
        block_size: must be smaller than N 

        return: list of coeffs list in which each of thems ranges from 0 to q - 1
        '''
        
        m = np.pad(m, pad_width=(0, (-len(m) % block_size)), mode='constant', constant_values=2) #very insecure padding 
        c = []
        for i in tqdm(range(0, len(m)//block_size)):
            c.append(self.encrypt_block(Poly(self.Rp, m[i * block_size: (i + 1) * block_size])).coeffs)
        
        return c
    
    def decrypt(self, enc: list, block_size = 96):
        '''
        enc: list of coeffs list in which each of thems ranges from 0 to q - 1
        block_size: must be smaller than N 

        return list of bit
        '''
        pad = lambda m : np.pad(m, pad_width=(0, (-len(m) % block_size)), mode='constant', constant_values=0)
        m = []
        for e in tqdm(enc):
            m.extend(pad(self.decrypt_block(Poly(self.Rq, e)).coeffs))
        
        print(m)
        i = 1
        while m[-i] == 2:
            i += 1
         
        return m[:-(i -  1)]
    
def test():

    #Moderate security parameter 
    N = 107
    p = 3
    q = 128
    df = 14
    dg = 12
    d = 5

    # rand = lambda : random.choice([1,2])
    # m = np.array([rand() + 1j*rand() for i in range(23)])
    # print(m)

   
    # ntru = NTRU(N, p, q, df, dg, d)
    # ntru.key_gen()
    # e = ntru.encrypt(m)
    # m_ = ntru.decrypt(e)
    # print('[+] Ciphertext:\n',e)
    # print('[+] Decrypted message:\n', m_ := ntru.decrypt(e))
  
    
    # if all([a == b for a, b in zip(m,m_)]):
    #     pass
    # else:
    
    #     print(m)
    #     print(m_)
    #     print('Wrong!')
    #     exit(0)

    image_path = 'lena.png'
    image = Image.open(image_path)

    # Convert the image to bytes
    with open(image_path, "rb") as image_file:
        image_bytes = image_file.read()

    # Convert bytes to binary
    m = list(map(int, bin(int(binascii.hexlify(image_bytes), 16))[2:]))[:200]
    
    ntru = NTRU(N, p, q, df, dg, d)
    ntru.key_gen()
    c = ntru.encrypt(m)
    m_ = ntru.decrypt(c)

    print(m)
    print(m_)
    assert(m == m_)

t = time.time()
for i in tqdm(range(1)):
    test()
print(time.time() - t)

