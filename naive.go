package fft

import (
  "math"
  "cmath"
  "fmt"
)
func init() {
  fmt.Printf("")
}
const J = complex(0.0, 1.0)

// Performs a DFT or a partial DFT on in, storing the output in out
// if start == 0 && stride == 1 the DFT is done on the entire array,
// otherwise it is done on every stride elements, starting at start.
func DFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  if N == 1 {
    out[start] = in[start]
    return
  }
  factor := -2 * math.Pi * complex(0,1) / complex(float64(N),0)
  for k := 0; k < N; k++ {
    out[start + k*stride] = 0
    for n := start; n < len(in); n += stride {
      out[start + k*stride] += in[n] * cmath.Exp(factor * complex(float64(k * (n / stride)), 0))
    }
  }
}

var n_factors []complex128
var kn_factors [][]complex128
const num_kn_factors = 64
const num_n_factors = num_kn_factors * num_kn_factors
func init() {
  n_factors = make([]complex128, num_n_factors)
  for i := range n_factors {
    n_factors[i] = -2 * math.Pi * J * complex(1.0 / float64(i), 0.0)
  }

  base_kn := make([]complex128, num_n_factors)
  kn_factors = make([][]complex128, num_kn_factors)
  for i := range kn_factors {
    kn_factors[i] = base_kn[i * num_kn_factors : (i+1) * num_kn_factors]
  }
  for i := range kn_factors {
    for j := range kn_factors {
      kfactor := n_factors[j] * complex(float64(i), 0.0)
      kn_factors[i][j] = cmath.Exp(kfactor)
    }
  }
}

func twiddle(f, N int) complex128 {
  d := -2 * math.Pi * J * complex(float64(f), 0) / complex(float64(N),0)
  m := cmath.Exp(d)
  return m
}

func butterfly(in,out []complex128, q, start, stride int) {
  // len(v) == p*q
  // TODO: this needs to be able to stride
  // start = 1, stride = 2
  //  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
  //     *     *     *     *     *     *     *     *     *     *
  // q = 5
  // p = 2
  // start + stride * (i + j*p)
  p := len(in) / q / stride
  for i := 0; i < p; i++ {
    for j := 0; j < q; j++ {
      out[start + stride * (i + j*p)] = 0
      for k := 0; k < q; k++ {
        out[start + stride * (i + j*p)] += in[start + stride * (i*q + k)] * twiddle((k * (i + j*p)), len(in) / stride)
      }
    }
  }
}

func factor(n int) []int {
  var f []int
  for i := 2; i*i <= n; i++ {
    for n%i == 0 {
      f = append(f, i)
      n /= i
      if n == 1 { return f }
    }
  }
  if n > 1 {
    f = append(f, n)
  }
  return f
}

func fftHelper(in,out,temp []complex128, factors []int, start,stride int) {
  if len(factors) == 1 {
    DFT(in, out, start, stride)
    return
  }
  factor := factors[0]
  factors = factors[1:]
  for i := 0; i < factor; i++ {
//    sft_helper(in, out, factors, stride*factor)
    fftHelper(in, out, temp, factors, start + stride*i, stride*factor)
  }
  copy(temp, out)

  butterfly(temp, out, factor, start, stride)
}

func FFT(in,out []complex128) {
  temp := make([]complex128, len(in))
  fftHelper(in, out, temp, factor(len(in)), 0, 1)
}
