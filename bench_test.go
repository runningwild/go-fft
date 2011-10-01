package fft_test

import (
  "testing"
  "fftw"
  "fft"
)

var N int = 1024*64

func BenchmarkFFTW(b *testing.B) {
  in := fftw.Alloc1d(N)
  out := fftw.Alloc1d(N)
  plan := fftw.PlanDft1d(in, out, fftw.Forward, fftw.Estimate)
  for i := 0; i < b.N; i++ {
    plan.Execute()
  }
}

func BenchmarkFFT(b *testing.B) {
  in := make([]complex128, N)
  out := make([]complex128, N)
  for i := 0; i < b.N; i++ {
    fft.FFT(in, out, 0, 1)
  }
}
