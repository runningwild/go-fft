package main

import (
  "fft"
  "fftw"
  "fmt"
  "math"
  "cmath"
)

func show(a []complex128) {
  for i,v := range a {
    re := real(v)
    im := real(v)
    if re < 1e-9 && re > -1e-9 { re = 0 }
    if im < 1e-9 && im > -1e-9 { im = 0 }
    fmt.Printf("%d:\t%2.4f\t%2.4f\n", i, re, im)
  }
  fmt.Printf("\n")
}

func twiddle(f, N int) complex128 {
  d := -2*math.Pi*fft.J * complex(float64(f), 0) / complex(float64(N),0)
  m := cmath.Exp(d)
  return m
}

func main() {
  N := 16
  in := make([]complex128, N)
  for i := range in {
    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * math.Pi * 2), 0)
  }
  out_fft := make([]complex128, N)
  out_fftw := make([]complex128, N)
  fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
  show(in)
  show(out_fftw)
  fft.DFT(in, out_fft, 0, 2)
  fft.DFT(in, out_fft, 1, 2)

  for i := range in {
    in[i] = complex(0,0)
  }
  no2 := len(in) / 2
  for i := 0; i < no2; i++ {
    in[2*i + 0] += out_fft[2*i + 0]
    in[2*i + 0] += out_fft[2*i + 1] * twiddle(i, len(in))
    in[2*i + 1] += out_fft[2*i + 0]
    in[2*i + 1] += out_fft[2*i + 1] * twiddle(i + no2, len(in))
  }
  show(in)
  fft.DFT(in[0:8], out_fft[0:8], 0, 1)
  fft.DFT(in[8: ], out_fft[8: ], 0, 1)
  show(out_fft)
}
