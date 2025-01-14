# InputParameters

To simplify and unify the creation of all simulation objects in MOOSE, all input parameters must be declared and populated through
a single "InputParameters" object. This ensures that every constructor in MOOSE is uniform and ensures that every object
can be created through MOOSE's Factory pattern. The InputParameters object is a collection of parameters, each one with
seperate attributes that can be used to finely control the behavior of the underlying object. For example, parameters can be
marked as required or optional, be provided with a default or not, and be used to enhance GUI interfaces that may be used to
programmatically generate input files for MOOSE.

The complete list of attributes for each input parameter:

!listing framework/include/utils/InputParameters.h start=struct Metadata end=Metadata & at

## Applying or Transferring Common Parameters

When building a custom Action, it is often useful to read in several parameters that will be used to directly set
parameters on objects being built by the custom Action. The `InputParameters` object contains a few useful methods for
applying or transferring common parameters to avoid several manual lines for setting these parameters. See the utility methods
and corresponding documentation here:

!listing framework/include/utils/InputParameters.h  start=BEGIN APPLY PARAMETER METHODS end=END APPLY PARAMETER METHODS