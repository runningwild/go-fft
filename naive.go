package fft

import (
  "math"
  "fmt"
  "cmath"
)

const j = complex(0.0, 1.0)

// Performs a DFT or a partial DFT on in, storing the output in out
// if start == 0 && stride == 1 the DFT is done on the entire array,
// otherwise it is done on every stride elements, starting at start.
func DFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  factor := -2 * math.Pi * complex(0,1) / complex(float64(N),0)
  for k := 0; k < N; k++ {
    out[k] = 0
    for n := start; n < len(in); n += stride {
      out[k] += in[n] * cmath.Exp(factor * complex(float64(k * (n / stride)), 0))
    }
  }
}

func FFT(in,out []complex128, start,stride int) {
  N := len(out)
  if N % 2 != 0 {
    DFT(in, out, start, stride)
    return
  }
  DFT(in, out[ : N/2], start, stride*2)
  DFT(in, out[N/2 : ], start+stride, stride*2)
  for i := range out { 
    fmt.Printf("%d: %2.2f\n", i, out[i])
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
