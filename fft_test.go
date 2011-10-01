package fft_test

import (
  . "gospec"
  "gospec"
  "math"
  "cmath"
  "fft"
  "fmt"
)

func NaiveSpec(c gospec.Context) {
  in := make([]complex128, 256*4096)
  for i := range in {
    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * 2 * math.Pi), 0)
  }
  fft.FFT(in)
  return
  c.Specify("Test basic DFT", func() {
    in := make([]complex128, 16)
    for i := range in {
      in[i] = complex(math.Cos(float64(i) / float64(len(in)) * 2 * math.Pi), 0)
    }
    out := make([]complex128, len(in))
    fft.DFT(in, out, 0, 1)
    for i := range out {
      mag,_ := cmath.Polar(out[i])
      if i == 1 || i == len(out) - 1 {
        c.Expect(mag, IsWithin(1e-9), float64(len(out))/2)
      } else {
        c.Expect(mag, IsWithin(1e-9), 0.0)
      }
    }
    out2 := fft.FFT(in)
    for i := range in {
      fmt.Printf("%d: %2.2f %2.2f\n", i, out[i], out2[i])
    }
  })
}