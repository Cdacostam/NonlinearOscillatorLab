(* Execute in a fresh WolframKernel with -noinit and -pacletbase. *)
archive = Environment["NOL_PACLET_PATH"];
base = Environment["NOL_PACLET_BASE"];
If[!StringQ[archive] || !FileExistsQ[archive] || !StringQ[base] ||
    ExpandFileName[base] =!= ExpandFileName[$UserBasePacletsDirectory],
  Print["Set NOL_PACLET_PATH and start the kernel with -pacletbase NOL_PACLET_BASE"]; Exit[2]
];
context = "NonlinearOscillatorLab`";
If[MemberQ[$Packages, context] || PacletFind["NonlinearOscillatorLab"] =!= {},
  Print["The verification requires a fresh paclet base"]; Exit[2]
];
installed = Check[PacletInstall[ExpandFileName[archive]], $Failed];
If[Head[installed] =!= PacletObject, Print["PacletInstall failed"]; Exit[1]];
If[Check[Needs[context]; MemberQ[$Packages, context], False] =!= True,
  Print["Needs failed to load the package"]; Exit[1]
];

(* Evaluate after Needs, so short symbol names resolve in the paclet context. *)
actual = ToExpression[
  "NOLPeriod[NOLModel[4, 0, 0], 2/5, 0, WorkingPrecision -> 40]"
];
If[ToExpression["NOLVersion"] =!= "0.2.1" || !NumberQ[actual] ||
    Abs[actual - N[Pi, 40]] > 10^-30,
  Print["Paclet smoke test failed: ", actual]; Exit[1]
];
Print["PASS: ", $Version, " | ", installed["Location"], " | harmonic period = ", actual];
