import numpy as np
from Zmod import Zmod
from Crypto.Util.number import isPrime
from sympy import mod_inverse

cut_off = 35 #this value can be modified to fasten multiplication 

class Poly:
    def __init__(self, n, coeffs: np.array = [0]):
        if isinstance(coeffs[0],int) or isinstance(coeffs[0], np.int64):
            self.coeffs = np.array([c % n for c in coeffs], dtype= np.int64)
        else: self.coeffs = coeffs
        self.n = np.int64(n)

    def trim(self):
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs = self.coeffs[:-1]
        return self
    
    def change_ring(self, n):
        self.n = np.int64(n)
        self.coeffs = np.array([c % n for c in self.coeffs], dtype=np.int64)

        return self
    
    def negate_inplace(self):
        for i in range(len(self.coeffs)):
            self.coeffs[i] = (-self.coeffs[i]) % self.n
        return self
    
    def negate(self):
        return Poly(self.n, self.coeffs.copy()).negate_inplace()
    
    def add_inplace(self, other):
        if isinstance(other, int) or isinstance(other, np.int64):
            self.coeffs[0] = (self.coeffs[0] + (other % self.n)) % self.n  # Convert int to ring element
            return self.trim()

        if not self.n == other.n:
            raise ValueError("Polynomials must be from the same ring.")

        max_len = max(len(self.coeffs), len(other.coeffs))
        new_coeffs = np.pad(self.coeffs, (0, max_len - len(self.coeffs)), mode='constant', constant_values=0)

        for i in range(len(other.coeffs)):
            new_coeffs[i] = (new_coeffs[i] + other.coeffs[i]) % self.n  # Modular addition

        self.coeffs = np.array(new_coeffs, dtype=np.int64)
        return self.trim()
    
    def __add__(self, other):
        return Poly(self.n, self.coeffs.copy()).add_inplace(other)
    
    def __sub__(self, other):
        return Poly(self.n, self.coeffs.copy()).add_inplace(other.negate())
    

    def mul_mod_convolution_ring(self, other, N):
        """Optimized Karatsuba multiplication for modular polynomials."""
        
        # Handle scalar multiplication
        if isinstance(other, (int, np.int64)):
            new_coeffs = (self.coeffs * other) % self.n
            return Poly(self.n, new_coeffs)

        if self.n != other.n:
            raise ValueError("Polynomials must be from the same ring.")

        def next_power_of_2(n):
            """Finds the next power of 2 greater than or equal to n."""
            return 1 if n == 0 else 2**(n - 1).bit_length()

        def karatsuba(a, b, N, mod):
            """Recursive Karatsuba polynomial multiplication modulo X^N - 1."""
            max_len = max(len(a), len(b))
            pow2_len = next_power_of_2(max_len)  # Ensure power-of-2 size

            # Pad polynomials to next power of 2
            a = np.pad(a, (0, pow2_len - len(a)), constant_values=0)
            b = np.pad(b, (0, pow2_len - len(b)), constant_values=0)

            if pow2_len <= cut_off:  # Base case: Direct multiplication
                c = np.zeros(2 * pow2_len - 1, dtype=np.int64)
                for i in range(len(a)):
                    for j in range(len(b)):
                        c[i + j] = (c[i + j] + a[i] * b[j]) % mod
            else:
                mid = pow2_len // 2  

               
                low_a, high_a = a[:mid], a[mid:]
                low_b, high_b = b[:mid], b[mid:]

                # Recursive Karatsuba calls
                c1 = karatsuba(low_a, low_b, N, mod)
                c3 = karatsuba(high_a, high_b, N, mod)

                
                sum_a = (low_a + high_a) % mod
                sum_b = (low_b + high_b) % mod
                c2 = karatsuba(sum_a, sum_b, N, mod)

                # Compute final middle term
                c2 = (c2 - c1 - c3) % mod  

                
                c = np.zeros(2 * pow2_len - 1, dtype=np.int64)
                c[:len(c1)] = c1
                c[mid:mid + len(c2)] = (c[mid:mid + len(c2)] + c2) % mod
                c[2 * mid:2 * mid + len(c3)] = (c[2 * mid:2 * mid + len(c3)] + c3) % mod

            
            if N != 0:
                for k in range(len(c) - 1, N - 1, -1):
                    c[k - N] = (c[k - N] + c[k]) % mod
                    c[k] = 0  

            return c[:N] if N != 0 else c

        # Perform multiplication and return result
        return Poly(self.n, karatsuba(self.coeffs, other.coeffs, N, self.n)).trim()



    def __mul__(self, other):
        return Poly(self.n, self.coeffs.copy()).mul_mod_convolution_ring(other, 0)
    
    def __eq__(self, other):
        return isinstance(other, Poly) and self.n == other.n and self.coeffs == other.coeffs
    
    def divmod_convol_inplace(self, N):
        n = len(self.coeffs)
        for k in range(n - 1, N - 1, -1):
           
            self.coeffs[k - N] = (self.coeffs[k - N] + self.coeffs[k]) % self.n
            self.coeffs[k] = 0
           
        
        return self.trim()
    
    def divmod_convol(self, N):
        return Poly(self.n, self.coeffs.copy()).divmod_convol_inplace(N)

    # based on this https://web.archive.org/web/20121005122331/http://www.securityinnovation.com/uploads/Crypto/NTRUTech014.pdf
    def inv_mod_xN_prime_inplace(self, N):
        """ Compute the modular inverse of a polynomial modulo x^N - 1 """
        #n has to be prime

        if not isPrime(self.n): return None
        k = 0
        b = Poly(self.n, [1])  # Identity polynomial
        c = Poly(self.n)  # Empty polynomial
        
        # Initialize g = x^N - 1
        a = [0] * (N + 1)
        a[N] = 1  
        a[0] = -1  
        g = Poly(self.n, a)
        try:
            while len(self.coeffs) > 1:
                # Remove leading zeros
                while self.coeffs[0] == 0 and len(self.coeffs) > 1:
                    self.coeffs = self.coeffs[1:]
                    c.coeffs = np.insert(c.coeffs, 0, 0)               
                    k += 1

                if len(self.coeffs) == 1:
                    if self.coeffs[0] == 0:
                        return None

                    f0_inv = np.int64(mod_inverse(self.coeffs[0], self.n))
                    b = b * f0_inv
                    shift = (N - k) % N
                    if shift > 0:
                        b.coeffs = np.insert(b.coeffs, 0, [0] * shift)
                    return b.divmod_convol_inplace(N)
                
                if len(self.coeffs) < len(g.coeffs):
                    self, g = g, self
                    b, c = c, b  

                u = (self.coeffs[0] * np.int64(mod_inverse(g.coeffs[0], self.n))) % self.n               
                self -= g * u
                b -= c * u
               
            
        except Exception as e:
                raise ValueError("No modular inverse exists") from e

    def inv_mod_xN_prime(self, N):
        f_ = Poly(self.n, self.coeffs.copy())
        return f_.inv_mod_xN_prime_inplace(N)
    
    def inv_mod_xN_prime_pow(self, N):
        # power of 2 is usually used in most cases
        # using hansel lemma to caculate
        q = 2 
        while self.n % q != 0: q += 1
        
        temp = self.n
        while temp != 1:
            if temp % q != 0: return None
            temp //= q

        a = self.coeffs
        b = Poly(q, a).inv_mod_xN_prime(N)
        
        while q < self.n:
            q = q * q
            b.change_ring(q)
            b = (Poly(q, a).mul_mod_convolution_ring(b, N).negate_inplace() + 2).mul_mod_convolution_ring(b,N)
            
        return b.change_ring(self.n)
        
    def center_lift(self):
        bound = self.n // 2

        new_coeffs = []
        for i in range(len(self.coeffs)):
            if self.coeffs[i] > bound: new_coeffs.append(self.coeffs[i] - self.n)
            else: new_coeffs.append(self.coeffs[i])
        
        return new_coeffs
    
    def __repr__(self):
        terms = []
        for d, coeff in enumerate(self.coeffs):
            if coeff != 0:  # Skip zero coefficients
                term = f"{coeff}x^{d}" if d > 0 else f"{coeff}"
                terms.append(term)

        poly_str = " + ".join(terms[::-1]) if terms else "0"
        return f"Polynomial {poly_str} in Ring of integers modulo {self.n}"



 
# a = Poly(7, [1, 1, 0, 0, 0, 1])
# b = Poly(7, [3, 1, 1, 0, 3, 0, 1, 0, 0, 1, 3])
# print(b.divmod_convol(3))
# print(a * np.int64(2))
# b_1 = b.inv_mod_xN_prime(4)
# # print(b_1)
# print(b + 2)

# f = Poly(32, [-1, 1, 1, 0, 1, 0, -1, 0, 0, 1, -1])

# print(f.inv_mod_xN_prime_pow(11))
# print(((f.inv_mod_xN_prime_pow(11)) * f).divmod_convol(11))


