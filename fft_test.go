package fft_test

import (
  . "gospec"
  "gospec"
  "math"
  "math/cmplx"
  "fft"
  "fftw"
  "fmt"
)

func FFTWSpec(c gospec.Context) {
  c.Specify("Check agains fftw", func() {
    N := 2*3*5*7*11 
    in := make([]complex128, N)
    out_fft := make([]complex128, N)
    out_fftw := make([]complex128, N)
    for i := range in {
      in[i] = complex(float64(i / 1000 - 100), float64(i) / 10)
    }

    fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
    fft.FFT(in, out_fft)

    tolerance := 1e-5
    for i := range out_fft {
      c.Expect(real(out_fft[i]), IsWithin(tolerance), real(out_fftw[i]))
      c.Expect(imag(out_fft[i]), IsWithin(tolerance), imag(out_fftw[i]))
    }
  })
}

func NaiveSpec(c gospec.Context) {
  in := make([]complex128, 8)
  for i := range in {
    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * 2 * math.Pi), 0)
  }
  verify := func(out []complex128, start,stride int, name string) {
    c.Specify(name, func() {
      for i := start; i < len(out); i += stride {
        mag,ang := cmplx.Polar(out[i])
        if i == start + stride || i == len(out) + start - stride {
          c.Expect(mag, IsWithin(1e-9), float64(len(out) / stride)/2)
          if real(out[i]) < 0 {
            if ang < 0 {
              c.Expect(ang, IsWithin(1e-9), -math.Pi)
            } else {
              c.Expect(ang, IsWithin(1e-9), math.Pi)
            }
          } else {
            c.Expect(ang, IsWithin(1e-9), 0.0)
          }
        } else {
          c.Expect(mag, IsWithin(1e-9), 0.0)
        }
      }
    })
  }

  c.Specify("Test basic DFT", func() {
    out := make([]complex128, len(in))

    fft.DFT(in, out, 0, 1)
    verify(out, 0, 1, "Start/Stride 0/1")

    in2 := make([]complex128, 2*len(in))
    out2 := make([]complex128, 2*len(in))
    for i := range in2 {
      in2[i] = in[i/2]
    }
    fft.DFT(in2, out2, 0, 2)
    verify(out2, 0, 2, "Start/Stride 0/2")


    fft.DFT(in2, out2, 1, 2)
    verify(out2, 1, 2, "Start/Stride 1/2")

    in5 := make([]complex128, 5*len(in))
    out5 := make([]complex128, 5*len(in))
    for i := range in5 {
      in5[i] = in[i/5]
    }
    for i := 0; i < 5; i++ {
      fft.DFT(in5, out5, i, 5)
      verify(out5, i, 5, fmt.Sprintf("Start/Stride %d/%d", i, len(in5)/len(in)))
    }
  })
}