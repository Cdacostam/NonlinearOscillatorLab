(* ::Package:: *)
(* NonlinearOscillatorLab 0.2.0 -- patch for the Wolfram Language package.

   Degenerate equilibria (V''(x) = 0, e.g. the inflection points x = +-1 of a = 1, b = -2, c = 1,
   where b^2 = 4 a c) are approached asymptotically, so their energy level is critical (infinite
   period), exactly as for a barrier. Version 0.1.0 only treated strict maxima (V'' < 0) as
   barriers. This patch adds the key "Degenerate" and includes these levels in
   "OnSeparatrixEnergy", matching the Python package v0.2.0.

   Usage:  Get[".../mathematica/NonlinearOscillatorLab.wl"];
           Get[".../mathematica/patches/degenerate_equilibria.wl"];

   STATUS: verified (three checks True with Mathematica 15.0.0) and merged into
   NonlinearOscillatorLab.wl in version 0.2.1; kept for the record. *)

Begin["NonlinearOscillatorLab`Private`"];

NOLClassifyEnergy[m_NOLModel, x0_, v0_, OptionsPattern[]] :=
 Module[{tol = OptionValue["Tolerance"], en, eq, mx, mn, dg},
  en = NOLEnergy[m, x0, v0];
  eq = NOLStationaryPoints[m];
  mx = {#, NOLPotential[m, #]} & /@ Select[eq, N[NOLCurvature[m, #]] < -tol &];
  mn = {#, NOLPotential[m, #]} & /@ Select[eq, N[NOLCurvature[m, #]] > tol &];
  dg = {#, NOLPotential[m, #]} & /@ Select[eq, Abs[N[NOLCurvature[m, #]]] <= tol &];
  <|"Energy" -> en, "Equilibria" -> eq, "Minima" -> mn, "Maxima" -> mx, "Degenerate" -> dg,
    "OnSeparatrixEnergy" -> AnyTrue[Join[mx, dg], Abs[N[en - #[[2]]]] <= tol (1 + Abs[N[#[[2]]]]) &]|>];

End[];

(* Checks (expected: True, True, True):
   FailureQ[Quiet@NOLPeriod[NOLModel[1, -2, 1], 0, 1/Sqrt[3]]]
   NOLClassifyEnergy[NOLModel[1, -2, 1], 0, 1/Sqrt[3]]["Degenerate"] === {{-1, 1/6}, {1, 1/6}}
   Abs[NOLPeriod[NOLModel[1, -2, 1], 1/2, 0] - 7.746275937880837] < 10^-12          (Python value)
*)
