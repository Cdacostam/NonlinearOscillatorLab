(* Run with NOL_PACLET_OUTPUT set to the archive output directory. *)
out = Environment["NOL_PACLET_OUTPUT"];
If[!StringQ[out] || out === "", Print["Set NOL_PACLET_OUTPUT"]; Exit[2]];
out = ExpandFileName[out];
If[!DirectoryQ[out], CreateDirectory[out, CreateIntermediateDirectories -> True]];
If[!DirectoryQ[out], Exit[2]];

source = DirectoryName[$InputFileName];
stagingRoot = CreateDirectory[];
staging = CreateDirectory[FileNameJoin[{stagingRoot, "NonlinearOscillatorLab"}]];
kernel = CreateDirectory[FileNameJoin[{staging, "Kernel"}]];
staged = Check[
  CopyFile[FileNameJoin[{source, "PacletInfo.wl"}], FileNameJoin[{staging, "PacletInfo.wl"}]];
  CopyFile[FileNameJoin[{source, "NonlinearOscillatorLab.wl"}], FileNameJoin[{kernel, "NonlinearOscillatorLab.wl"}]];
  CopyFile[FileNameJoin[{DirectoryName[source], "LICENSE"}], FileNameJoin[{staging, "LICENSE"}]];
  True,
  False
];
archive = If[TrueQ[staged], Check[CreatePacletArchive[staging], $Failed], $Failed];
If[!StringQ[archive] || !FileExistsQ[archive],
  DeleteDirectory[stagingRoot, DeleteContents -> True]; Print["Paclet build failed"]; Exit[1]
];
target = FileNameJoin[{out, FileNameTake[archive]}];
If[FileExistsQ[target], Print["Archive already exists at ", target]; DeleteDirectory[stagingRoot, DeleteContents -> True]; Exit[2]];
result = Check[CopyFile[archive, target], $Failed];
DeleteDirectory[stagingRoot, DeleteContents -> True];
If[result === $Failed || !FileExistsQ[target], Exit[1]];
Print[target];
