# Saving `tck`

Once a tractogram has been computed, it can be exported to `tck` file using `nibabel`. A convenience function is provided to ease this operation:

```julia
import Tractography as TG

# ... model definition skipped
streamlines, tract_length = TG.sample(model, alg_tracking, seeds)

# save to tck file
using PythonCall
TG.save_streamlines("my_stream.tck", streamlines, tract_length)
```