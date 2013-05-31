plantecophys package notes
==========================


## To do list

### High priority

* T response parameters should be a list?
 * photosyn() (vectorized) calculates parameters based on T response
 * sends it to photosynF() (not vectorized)
 * avoids calculating it many times!
* get MESO from GasExchangeR photosyn (BM's edits)
* `photosyn()` and `Aci()` will have substantial overlap, remove overlap as much as possible.


### Low priority

* Roxygen
* no Tuzet, yet
* no Penman-Monteith (as MAESTRA)
* more unit conversion utilities?


## Notes
* `photosyn2()` is faster than `GasExchangeR::photosyn()` (see doc in subdir)


