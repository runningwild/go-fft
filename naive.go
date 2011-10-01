package fft

import (
  "math"
  "cmath"
)

const j = complex(0.0, 1.0)

// Performs a DFT or a partial DFT on in, storing the output in out
// if start == 0 && stride == 1 the DFT is done on the entire array,
// otherwise it is done on every stride elements, starting at start.
func DFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  factor := -2 * math.Pi * complex(0,1) / complex(float64(N),0)
  for i := range out {
    out[i] = 0
  }
  index := 0
  for k := range out {
    for n := start; n < len(in); n += stride {
      out[k] += in[n] * cmath.Exp(factor * complex(float64((k * n) / stride), 0))
    }
    index++
  }
}

func FFT(in []complex128) (out []complex128) {
  out = make([]complex128, len(in))
  N := len(out)
  if N % 2 != 0 {
    DFT(in, out, 0, 1)
    return
  }
  i0 := make([]complex128, N/2)
  i1 := make([]complex128, N/2)
  for i := range i0 {
    i0[i] = in[2*i]
    i1[i] = in[2*i+1]
  }
  o1 := FFT(i0)
  o2 := FFT(i1)
  for i := range o1 {
    out[i] = o1[i]
    out[i+N/2] = o2[i]
  }
  factor := -2 * math.Pi * j * complex(1.0 / float64(N), 0.0)
  for k := 0; k < N/2; k ++ {
    t := out[k]
    kfactor := factor * complex(float64(k), 0.0)
//    kmnfactor := factor * complex(float64(k-N/2), 0.0)
    out[k]     = t + cmath.Exp(kfactor) * out[k + N/2]
    out[k+N/2] = t - cmath.Exp(kfactor) * out[k + N/2]
  }
  return
}
