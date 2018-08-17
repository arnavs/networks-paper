# networks.jl

Provides an API around the model described in Fogli and Veldkamp, "Germs, Social Networks, and Growth." (forthcoming).

date: 2017-09-24

## The Model

This is an endogenous model of social network formation. It's "endogenous" because, while the social graph influences the economic state (e.g., governs how technology and disease spread from person to person), these things also influence the evolution of the social graph. The alternative (an exogenous model) would be to peg one of those two to a fixed process (which could be stochastic).

For an example of the social graphs we're talking about, see: 

![Alt text](output-CanonicalModel.png?raw=true "Canonical Model Output")

Note: This paper might not be forthcoming anymore. The authors decided to stick with MATLAB, which means I didn't have much continuing involvement. If you're interested, you can reach out to Fogli or Veldkamp. 

## This Repo 

The meat of the model is in `networksBackend.jl`. There's code to inspect the results and grab relevant data in `...Reporting.jl`, code to run a simulation or two in `...Simulation.jl`, and then the output above from one set of parameters in the png file. 

## Contact 

If you have an interest in this Julia code specifically, reach out to me at [arnav.sood@ubc.ca](mailto:arnav.sood@ubc.ca).
