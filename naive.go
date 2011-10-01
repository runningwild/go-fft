package fft

import (
  "math"
  "cmath"
)

const J = complex(0.0, 1.0)

// Performs a DFT or a partial DFT on in, storing the output in out
// if start == 0 && stride == 1 the DFT is done on the entire array,
// otherwise it is done on every stride elements, starting at start.
func DFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  if N == 1 {
    out[0] = in[start]
    return
  }
  factor := -2 * math.Pi * complex(0,1) / complex(float64(N),0)
  for k := 0; k < N; k++ {
    out[k] = 0
    for n := start; n < len(in); n += stride {
      out[k] += in[n] * cmath.Exp(factor * complex(float64(k * (n / stride)), 0))
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

func FFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  if N == 1 {
    out[0] = in[start]
    return
  }
  if N % 2 != 0 {
    DFT(in, out, start, stride)
    return
  }
  FFT(in, out[ : N/2], start, stride*2)
  FFT(in, out[N/2 : ], start+stride, stride*2)

  var factor complex128
  if N < num_n_factors {
    factor = n_factors[N]
  } else {
    factor = -2 * math.Pi * J * complex(1.0 / float64(N), 0.0)
  }

  var knfactor complex128
  for k := 0; k < N/2; k ++ {
    t := out[k]

    if N < num_kn_factors {
      knfactor = kn_factors[k][N]
    } else {
      knfactor = cmath.Exp(factor * complex(float64(k), 0.0))
    }
    term := knfactor * out[k + N/2]
    out[k]     = t + term
    out[k+N/2] = t - term
  }
  return
}
