#!/usr/bin/env wolframscript
(* T1: three checks of patches/degenerate_equilibria.wl.
   A) ORIGINAL v0.2.0 package + patch;  B) merged v0.2.1 package (no patch).
   Each case runs in its own fresh kernel; the checks are parsed only AFTER the package is loaded. *)
$dir = DirectoryName[$InputFileName];
$checks = "{FailureQ[Quiet@NOLPeriod[NOLModel[1, -2, 1], 0, 1/Sqrt[3]]],
  NOLClassifyEnergy[NOLModel[1, -2, 1], 0, 1/Sqrt[3]][\"Degenerate\"] === {{-1, 1/6}, {1, 1/6}},
  Abs[NOLPeriod[NOLModel[1, -2, 1], 1/2, 0] - 7.746275937880837] < 10^-12}";
caseCode[files_] := StringJoin[
  "Get /@ ", ToString[files, InputForm], ";\n",
  "res = Check[", $checks, ", \"MESSAGES\"];\n",
  "Print[\"checks: \", ToString[res, InputForm]];\n",
  "Print[\"NOLPeriod[NOLModel[1,-2,1], 1/2, 0] = \", ToString[NOLPeriod[NOLModel[1, -2, 1], 1/2, 0], InputForm]];\n",
  "Print[\"Degenerate = \", ToString[NOLClassifyEnergy[NOLModel[1, -2, 1], 0, 1/Sqrt[3]][\"Degenerate\"], InputForm]];\n",
  "Print[\"NOLVersion = \", NOLVersion];\n"];
runCase[files_] := Module[{f = FileNameJoin[{$TemporaryDirectory, "nol_t1_" <> ToString[RandomInteger[10^9]] <> ".wls"}], r},
  Export[f, caseCode[files], "Text"];
  r = RunProcess[{First[$CommandLine], "-script", f}];
  DeleteFile[f];
  StringTrim[r["StandardOutput"] <> r["StandardError"]]];
txt = StringJoin[
  "$Version: ", $Version, "\n\n",
  "A) original v0.2.0 + patches/degenerate_equilibria.wl\n",
  runCase[{FileNameJoin[{$dir, "original_v0.2.0", "NonlinearOscillatorLab.wl"}],
           FileNameJoin[{$dir, "original_v0.2.0", "patches", "degenerate_equilibria.wl"}]}], "\n\n",
  "B) merged v0.2.1 (no patch)\n",
  runCase[{FileNameJoin[{$dir, "v0.2.1", "NonlinearOscillatorLab.wl"}]}], "\n"];
Export[FileNameJoin[{$dir, "T1_checks_output.txt"}], txt, "Text"];
Print[txt];
