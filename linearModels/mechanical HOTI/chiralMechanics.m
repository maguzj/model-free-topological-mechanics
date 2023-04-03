(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Jan 6, 2022 *)

BeginPackage["chiralMechanics`"]
(* Exported symbols added here with SymbolName::usage *) 

eqSystem::usage = "Gives a graph in which the nodes are beads and the edges are springs."
eqSystem::badarg = "The functions expects three arguments: keysBeads_List,pairKeysBeads_List,dsPositions_Association"

pertSystem::usage = "Gives a chiral graph of displacementes and elongations."
pertSystem::badarg = "The function expects four arguments: graphBeads_Graph, dsBeadIdentity_Association, dsPositions_Association, dsPivotPositions_Association."

pertSystem3D::usage = "Gives a chiral graph of displacements and elongations in 3D"
pertSystem3D::badarg = "wrong format/arguments"

dElongation::usage = "Gives the variation of elongation to first order due to displacements."
dElongation::badarg = "wrong format/arguments"

(*assignName::usage = "a function to assign names for the beads in a unit cell"
assignName::badarg = "wrong format/arguments"

it::usage = "iterator"
it::badarg = "wrong format/arguments"*)

tesselate::usage = "Returns {keysBeads, pairKeysBeads, identity, posBeads, posPivots} for a crystal made of {nx,ny,nz} unit cells"
tesselate::badarg = "wrong format/arguments"

beadAndSpringGraphic2D::usage = "Returns a graphics object showin rotors, beads, and springs. Input: output of tesselate"
beadAndSpringGraphic2D::badarg = "wrong format/arguments"

beadAndSpringGraphic3D::usage = "Returns a graphics object showin rotors, beads, and springs. Input: output of tesselate"
beadAndSpringGraphic3D::badarg = "wrong format/arguments"

blochCMatrix::usage = "returns the associated chiral Bloch hamiltonian. Input from the models package"
blochCMatrix::badarg = "wrong format/arguments"

(*highlightGraph::usage = "returns the perturbed graph with opacity given by the input state"
highlightGraph::badarg = "wrong format/arguments"*)

rectangle::usage = "give two middle points and wait"
rectangle::badarg = "ups"

Begin["`Private`"] (* Begin Private Context *) 

eqSystem[keysBeads_List,pairKeysBeads_List,dsPositions_Association] :=
	Graph[#[[1]] \[UndirectedEdge] #[[2]] & /@ pairKeysBeads, VertexCoordinates -> (# -> dsPositions[#]&/@keysBeads)]
eqSystem[___] := (Message[eqSystem::badarg];$Failed)

dElongation[r1_List,r2_List,\[Delta]r1_List,\[Delta]r2_List]/;(Equal @@ (Length[#] & /@ {r1, r2, \[Delta]r1, 
       \[Delta]r2})):=(r1-r2) . (\[Delta]r1-\[Delta]r2)/Norm[r1-r2]
dElongation[___] := (Message[dElongation::badarg];$Failed)

pertSystem[eqSystem_Graph, dsBeadIdentity_Association, dsPositions_Association, dsPivotPositions_Association] :=
    Module[ {weightedEdges, vertices, vertexStyle},
        {weightedEdges, vertices, vertexStyle} = 
        Transpose[
        Module[ {disp1, dofs1, labels1, disp2, labels2, dofs2, \[Delta]e, 
        pos1, pos2, \[Delta]elabel},
            Table[{{disp1, dofs1, labels1, pos1}, {disp2, dofs2, labels2, 
                pos2}} = 
              If[ dsBeadIdentity[edge[[#]]] == 
                  1,
                  {RotationMatrix[\[Pi]/2] . 
                  Normalize[(dsPositions[edge[[#]]] - 
                    dsPivotPositions[
                     edge[[#]]])] \[Delta]\[Theta][#], \
                  {\[Delta]\[Theta][#]}, {"\[Delta]\[Theta]" <> 
                  edge[[#]]}, {dsPositions[
                  edge[[#]]]}},
                  {{\[Delta]x[#], \[Delta]y[#]}, \
                  {\[Delta]x[#], \[Delta]y[#]}, {"\[Delta]x" <> edge[[#]], 
                  "\[Delta]y" <> edge[[#]]}, {dsPositions[edge[[#]]], 
                  dsPositions[edge[[#]]] + {0.1, 0.1}}}
              ] & /@ {1, 2};
                  \[Delta]e = (dElongation @@ {dsPositions[#1[[1]]], 
                        dsPositions[#1[[2]]], disp1, disp2}) &@edge;
                  \[Delta]elabel = edge[[1]] <> edge[[2]];
                  {Flatten@{MapThread[
                      Annotation[\[Delta]elabel \[UndirectedEdge] #1, 
                        EdgeWeight -> #2] &, {labels1, 
                       D[\[Delta]e, #] & /@ dofs1}], 
                     MapThread[
                      Annotation[\[Delta]elabel \[UndirectedEdge] #1, 
                        EdgeWeight -> #2] &, {labels2, 
                       D[\[Delta]e, #] & /@ dofs2}]}, 
                   Flatten@{\[Delta]elabel -> (dsPositions[edge[[1]]] + 
                         dsPositions[edge[[2]]])/2, Thread[labels1 -> pos1], 
                     Thread[labels2 -> pos2]}, 
                   Flatten@{\[Delta]elabel -> Blue, Thread[labels1 -> Red], 
                     Thread[labels2 -> Red]}}, {edge, EdgeList[eqSystem]}]
        ]];
        Graph[Flatten[weightedEdges], VertexCoordinates -> Flatten[vertices],
          VertexStyle -> Flatten[vertexStyle]]
    ]
pertSystem[___] := (Message[pertSystem::badarg];$Failed)


pertSystem3D[eqSystem_Graph, dsBeadIdentity_Association, 
  dsPositions_Association, dsPivotPositions_Association] :=
    Module[ {weightedEdges, vertices, 
      vertexStyle},
        {weightedEdges, vertices, vertexStyle} = 
        Transpose[
        Module[ {disp1, dofs1, labels1, disp2, labels2, dofs2, \[Delta]e, 
        pos1, pos2, \[Delta]elabel},
            Table[{{disp1, dofs1, labels1, pos1}, {disp2, dofs2, labels2, 
                pos2}} = If[
                 (*rotor case*) dsBeadIdentity[edge[[#]]] == 1,
                             {Total@
                               MapThread[#1 #2 &, \
                             {(Orthogonalize[{Normalize[(dsPositions[edge[[#]]] - 
                                       dsPivotPositions[edge[[#]]])]}~Join~
                  IdentityMatrix[3]][[{2, 
                                     3}]]), {\[Delta]\[Theta][#], \
\[Delta]\[CurlyPhi][#]}}],
                              {\[Delta]\[Theta][#], \[Delta]\[CurlyPhi][#]},
                              {"\[Delta]\[Theta]" <> edge[[#]], 
                               "\[Delta]\[CurlyPhi]" <> edge[[#]]},
                              {dsPositions[edge[[#]]], 
                               dsPositions[edge[[#]]] + {0.1, 0.1, 0.1}}},
                             (*free bead case*)
                             {{\[Delta]x[#], \[Delta]y[#], \
\[Delta]z[#]},
                              {\[Delta]x[#], \[Delta]y[#], \[Delta]z[#]},
                              {"\[Delta]x" <> edge[[#]], "\[Delta]y" <> edge[[#]], 
                               "\[Delta]z" <> edge[[#]]},
                              {dsPositions[edge[[#]]], 
                               dsPositions[edge[[#]]] + {0.1, 0.1, 0.1}, 
                               dsPositions[edge[[#]]] - {0.1, 0.1, 0.1}}}
                         ] & /@ {1, 2};
                  \[Delta]e = (dElongation @@ {dsPositions[#1[[1]]], 
                        dsPositions[#1[[2]]], disp1, disp2}) &@edge;
                  \[Delta]elabel = edge[[1]] <> edge[[2]];
                  {Flatten@{MapThread[
                      Annotation[\[Delta]elabel \[UndirectedEdge] #1, 
                        EdgeWeight -> #2] &, {labels1, 
                       D[\[Delta]e, #] & /@ dofs1}], 
                     MapThread[
                      Annotation[\[Delta]elabel \[UndirectedEdge] #1, 
                        EdgeWeight -> #2] &, {labels2, 
                       D[\[Delta]e, #] & /@ dofs2}]}, 
                   Flatten@{\[Delta]elabel -> (dsPositions[edge[[1]]] + 
                         dsPositions[edge[[2]]])/2, Thread[labels1 -> pos1], 
                     Thread[labels2 -> pos2]}, 
                   Flatten@{\[Delta]elabel -> Blue, Thread[labels1 -> Red], 
                     Thread[labels2 -> Red]}}, {edge, EdgeList[eqSystem]}]
        ]];
        Graph[Flatten[weightedEdges], 
         VertexCoordinates -> Flatten[vertices], 
         VertexStyle -> Flatten[vertexStyle]]
    ]
pertSystem3D[___] := (Message[pertSystem3D::badarg];$Failed)


(*assignName[s_String, i_String] := s <> i;
assignName[s_String, l_List] := assignName[s, StringRiffle[l, ","]];
assignName[s_String, i_Integer] := assignName[s, ToString[i]];
assignName[___] := (Message[assignName::badarg];$Failed)*)


Clear[it];
it[f_Symbol, l_List] /; Length@l == 1 := Table[f[nx], {nx, 1, l[[1]]}];
it[f_Symbol, l_List] /; Length@l == 2 := 
  Table[f @@ {nx, ny}, {nx, 1, l[[1]]}, {ny, 1, l[[2]]}];
it[f_Symbol, l_List] /; Length@l == 3 := 
  Table[f @@ {nx, ny, nz}, {nx, 1, l[[1]]}, {ny, 1, l[[2]]}, {nz, 1, 
    l[[3]]}];
it[___] := $Failed;


Clear[tesselate];
tesselate[nCells_List, keysBeadsUC_Symbol, pairKeysUC_Symbol, 
    identityUC_List, posBeadsUC_List, translationOp_Symbol, 
    posPivotsUC_List] :=
    Module[ {keysBeads, pairKeysBeads, posBeads, identity, posPivots},
    (**)
        keysBeads = Flatten[it[keysBeadsUC, nCells]];
        pairKeysBeads = Flatten[it[pairKeysUC, nCells], Length@nCells];
        (*delete extra cases that are not in keysBeads*)
        pairKeysBeads = 
         DeleteCases[pairKeysBeads, 
          x_ /; (!MemberQ[keysBeads, x[[2]]] \[Or] !MemberQ[keysBeads, x[[1]]])];
          
        keysBeads = ToString[#] & /@ keysBeads;
        pairKeysBeads = Map[ToString[#] &, pairKeysBeads, {2}];
        
        identity = Flatten@ConstantArray[identityUC, nCells];
        posBeads = Flatten[With[ {len = Length@posBeadsUC},
                               (posBeadsUC + ConstantArray[#, len]) & /@ 
                                Flatten[it[translationOp, nCells], Length@nCells]
                           ], 1];
        posPivots = Flatten[With[ {len = Length@posPivotsUC},
                                (posPivotsUC + ConstantArray[#, len]) & /@ 
                                 Flatten[it[translationOp, nCells], Length@nCells]
                            ], 1];
        {keysBeads, pairKeysBeads, identity, posBeads, posPivots}
    ];
tesselate[___] := (Message[tesselate::badarg];$Failed);


Clear[rectangle];
rectangle[{p1_List, p2_List}, width_Real /; width > 0] :=
    Module[ {t = (p2 - p1)/Norm[p2 - p1], n},
        n = RotationMatrix[\[Pi]/2] . t;
        Polygon[{p1 + width n, p1 - width n, p2 - width n, p2 + width n}]
    ]
rectangle[___] :=
    $Failed;

Clear[spring];
spring[{p1_List, p2_List}, width_Real /; width > 0, 
  nCoils_Integer : 10] :=
    Module[ {t = (p2 - p1)/Norm[p2 - p1], n},
        n = RotationMatrix[\[Pi]/2] . t;
        Line[p1 + #/nCoils (p2 - p1) + (-1)^# width n & /@ Range[0, nCoils]]
    ]
spring[___] :=
    $Failed;
    
Clear[spring3D];
spring3D[{p1_List, p2_List}, spiralRadius_Real /; spiralRadius > 0, 
  tubeRadius_ : 0.02, nSpirals_ : 5] := 
 Module[{t = (p2 - p1)/Norm[p2 - p1], n1, n2},
  {n1, n2} = 
   Select[Orthogonalize[{t}~Join~IdentityMatrix[3]], 
     Norm[#] > 0.5 &][[{2, 3}]];
  Tube[BSplineCurve[
    p1 + # (p2 - p1) + 
       spiralRadius (Sin[nSpirals 2 \[Pi] #] n1 + 
          Cos[nSpirals 2 \[Pi] #] n2) & /@ Range[0, 1, 0.05]], 
   tubeRadius]]
spring3D[___] := $Failed;


beadAndSpringGraphic2D[{keysBeads_List, pairKeysBeads_List, identity_List, posBeads_List, 
    posPivots_List}, 
   radius_Real : 
    0.1] /; (Equal @@ (Length[#] & /@ {keysBeads, identity, posBeads, 
       posPivots})) := Module[
  {diskList, posBeadsRotors, posPivotsRotors, rotorList, 
   dsPositionsBeads, springList},
  
  (*disks*)
  diskList = Disk[#, radius] & /@ posBeads;
  (*rotors*)
  posBeadsRotors = Pick[posBeads, identity, 1];
  posPivotsRotors = Pick[posPivots, identity, 1];
  rotorList = MapThread[
    rectangle[{#2, #1 - (radius) (#1 - #2)/Norm[#1 - #2]}, radius/
      2] &, {posBeadsRotors, posPivotsRotors}];
  (*springs*)
  
  dsPositionsBeads = AssociationThread[keysBeads, posBeads];
  springList = 
   With[{p1 = dsPositionsBeads[#[[1]]], p2 = dsPositionsBeads[#[[2]]]},
      (*rectangle[{p1+radius(p2-p1)/Norm[p2-p1],p2-radius(p2-p1)/Norm[
      p2-p1]},radius]*)
      spring[{p1, p2}, radius/2]
      ] & /@ pairKeysBeads;
  
  (*output*)
  
  Graphics[{Opacity[0.7], Black, rotorList, Blue, Thick, springList, 
    Red, diskList}]
  
  ]
beadAndSpringGraphic2D[___] := (Message[beadAndSpringGraphic2D::badarg];$Failed)

beadAndSpringGraphic3D[{keysBeads_List, pairKeysBeads_List, 
    identity_List, posBeads_List, posPivots_List}, 
   radius_Real : 
    0.1] /; (Equal @@ (Length[#] & /@ {keysBeads, identity, posBeads, 
       posPivots})) := 
 Module[{sphereList, posBeadsRotors, posPivotsRotors, rotorList, 
   dsPositionsBeads, springList},
  (*spheres*)sphereList = Sphere[#, radius] & /@ posBeads;
  (*rotors*)posBeadsRotors = Pick[posBeads, identity, 1];
  posPivotsRotors = Pick[posPivots, identity, 1];
  rotorList = 
   MapThread[
    Cylinder[{#2, #1 - (radius) (#1 - #2)/Norm[#1 - #2]}, 
      radius/2] &, {posBeadsRotors, posPivotsRotors}];
  (*springs*)
  
  dsPositionsBeads = AssociationThread[keysBeads, posBeads];
  springList = 
   With[{p1 = dsPositionsBeads[#[[1]]], 
       p2 = dsPositionsBeads[#[[2]]]},(*rectangle[{p1+radius(p2-p1)/
      Norm[p2-p1],p2-radius(p2-p1)/Norm[p2-p1]},radius]*)
      spring3D[{p1, p2}, radius/2]] & /@ pairKeysBeads;
  (*output*)
  Graphics3D[{Opacity[0.7], Black, rotorList, Blue, Thick, springList,
     Red, sphereList}]]
beadAndSpringGraphic3D[___] := (Message[beadAndSpringGraphic3D::badarg];$Failed)

blochCMatrix[
   kVector_List /; Length[kVector] == 2, {dimCrystal_, 
    dimEmbedding_Integer}, keysBeadsUC_, pairKeysUC_, identityUC_, 
   posBeadsUC_, translationOp_, 
   posPivotsUC_] /; (Dimensions[posBeadsUC][[2]] == 
     dimEmbedding \[And] Length@kVector <= dimEmbedding) := 
 Module[{dsPositions, dsPivotPositions, dsBeadIdentity, cMatrix, 
   disp1, disp2, i, dofs, fPos, fPiv, fId}, 
  fPos[arg__] := 
   AssociationThread[keysBeadsUC@arg, 
    Table[(translationOp@arg)[[1]] + pos, {pos, posBeadsUC}]];
  fPiv[arg__] := 
   AssociationThread[keysBeadsUC@arg, 
    Table[(translationOp@arg)[[1]] + pos, {pos, posPivotsUC}]];
  fId[arg__] := AssociationThread[keysBeadsUC@arg, identityUC];
  i = 1;
  dsPositions = Association@it[fPos, ConstantArray[i + 1, dimCrystal]];
  dsPivotPositions = 
   Association@it[fPiv, ConstantArray[i + 1, dimCrystal]];
  dsBeadIdentity = 
   Association@it[fId, ConstantArray[i + 1, dimCrystal]];
  (*dsPositions=Association@Table[AssociationThread[keysBeadsUC[ii],
  translationOp[ii][[1]]+#&/@posBeadsUC],{ii,{i,i+1}}];
  dsPivotPositions=Association@Table[AssociationThread[keysBeadsUC[
  ii],translationOp[ii][[1]]+#&/@posPivotsUC],{ii,{i,i+1}}];
  dsBeadIdentity=Association@Table[AssociationThread[keysBeadsUC[ii],
  identityUC],{ii,{i,i+1}}];*)(*degrees of freedom*)
  dofs = 
   Flatten[
    If[dsBeadIdentity[#] == 
        1, {\[Delta]\[Theta][Head[#]]}, {\[Delta]x[
         Head[#]], \[Delta]y[Head[#]]}] & /@ (keysBeadsUC @@ 
       ConstantArray[i, dimCrystal])];
  (*For each pair of of atoms determine the displacement times the \
bloch factor*)
  cMatrix = 
   Table[{{disp1}, {disp2}} = 
     If[dsBeadIdentity[edge[[#]]] == 
         1, {RotationMatrix[\[Pi]/2] . 
           Normalize[(dsPositions[edge[[#]]] - 
              dsPivotPositions[edge[[#]]])] \[Delta]\[Theta][
           Head[
            edge[[#]]]] (E^(-I translationOp @@ (Level[edge[[#]], 1] -
                    i) . kVector))[[1]]}, {{\[Delta]x[
            Head[edge[[#]]]], \[Delta]y[
            Head[
             edge[[#]]]]} (E^(-I translationOp @@ (Level[edge[[#]], 
                    1] - i) . kVector))[[1]]}] & /@ {1, 2};
    (*Compute the elongation using the displacements*)\[Delta]e = 
     Simplify@(dElongation @@ {dsPositions[#1[[1]]], 
           dsPositions[#1[[2]]], disp1, disp2}) &@edge;
    D[\[Delta]e, #] & /@ dofs, {edge, 
     pairKeysUC @@ ConstantArray[i, dimCrystal]}];
  (*h=Refine[
  ArrayFlatten[{{0,(cMatrix)\[HermitianConjugate]},{cMatrix,0}}],
  Element[kVector,Reals]]*)cMatrix]
blochCMatrix[
   kVector_List /; Length[kVector] == 3, {dimCrystal_, 
    dimEmbedding_Integer}, keysBeadsUC_, pairKeysUC_, identityUC_, 
   posBeadsUC_, translationOp_, 
   posPivotsUC_] /; (Dimensions[posBeadsUC][[2]] == 
     dimEmbedding \[And] Length@kVector <= dimEmbedding) := 
 Module[{dsPositions, dsPivotPositions, dsBeadIdentity, cMatrix, 
   disp1, disp2, i, dofs, fPos, fPiv, fId}, 
  fPos[arg__] := 
   AssociationThread[keysBeadsUC@arg, 
    Table[(translationOp@arg)[[1]] + pos, {pos, posBeadsUC}]];
  fPiv[arg__] := 
   AssociationThread[keysBeadsUC@arg, 
    Table[(translationOp@arg)[[1]] + pos, {pos, posPivotsUC}]];
  fId[arg__] := AssociationThread[keysBeadsUC@arg, identityUC];
  i = 1;
  dsPositions = Association@it[fPos, ConstantArray[i + 1, dimCrystal]];
  dsPivotPositions = 
   Association@it[fPiv, ConstantArray[i + 1, dimCrystal]];
  dsBeadIdentity = 
   Association@it[fId, ConstantArray[i + 1, dimCrystal]];
  (*dsPositions=Association@Table[AssociationThread[keysBeadsUC[ii],
  translationOp[ii][[1]]+#&/@posBeadsUC],{ii,{i,i+1}}];
  dsPivotPositions=Association@Table[AssociationThread[keysBeadsUC[
  ii],translationOp[ii][[1]]+#&/@posPivotsUC],{ii,{i,i+1}}];
  dsBeadIdentity=Association@Table[AssociationThread[keysBeadsUC[ii],
  identityUC],{ii,{i,i+1}}];*)(*degrees of freedom*)
  dofs = 
    Flatten[
     If[dsBeadIdentity[#] == 
         1, {\[Delta]\[Theta][Head[#]], \[Delta]\[CurlyPhi][
          Head[#]]}, {\[Delta]x[Head[#]], \[Delta]y[
          Head[#]], \[Delta]z[Head[#]]}] & /@ (keysBeadsUC @@ 
        ConstantArray[i, dimCrystal])];
  (*For each pair of of atoms determine the displacement times the \
bloch factor*)
  cMatrix = 
   Table[{{disp1}, {disp2}} = 
     If[dsBeadIdentity[edge[[#]]] == 1, 
        {Total[
            MapThread[#1 #2 &, {(Sort[
                 Orthogonalize[{Normalize[(dsPositions[edge[[#]]] - 
                    dsPivotPositions[edge[[#]]])]}~Join~
                    IdentityMatrix[3]][[2 ;;]]][[{2, 
                  3}]]), {\[Delta]\[Theta][
                Head[edge[[#]]]], \[Delta]\[CurlyPhi][
                Head[
                 edge[[#]]]]}}]] (E^(-I translationOp @@ (Level[
                    edge[[#]], 1] - i) . kVector))[[1]]}, {{\[Delta]x[
            Head[edge[[#]]]], \[Delta]y[Head[edge[[#]]]], \[Delta]z[
            Head[
             edge[[#]]]]} (E^(-I translationOp @@ (Level[edge[[#]], 
                    1] - i) . kVector))[[1]]}] & /@ {1, 2};
    (*Compute the elongation using the displacements*)
    \[Delta]e = 
      Simplify@(dElongation @@ {dsPositions[#1[[1]]], 
            dsPositions[#1[[2]]], disp1, disp2}) &@edge;
    D[\[Delta]e, #] & /@ dofs, {edge, 
     pairKeysUC @@ ConstantArray[i, dimCrystal]}];
  (*h=Refine[
  ArrayFlatten[{{0,(cMatrix)\[HermitianConjugate]},{cMatrix,0}}],
  Element[kVector,Reals]]*)cMatrix]
blochCMatrix[___] := (Message[blochCMatrix::badarg];$Failed)



(*this implementation follows the fact that the adjancencyMatrix has \
the same order of columns ands rows as given by VertexList. See \
Documentation*)
(*highlightGraph[graph_Graph,state_List,baseOpacity_:0.2]:=
Graph[graph, 
 VertexStyle -> 
  MapThread[#1 -> 
     Directive[EdgeForm[None], 
      If[StringPart[#1, 1] == "\[Delta]", Red, Blue], 
      Opacity[baseOpacity + #2]] &, {VertexList[graph], 
    state}]]
highlightGraph[___] := (Message[highlightGraph::badarg];$Failed)*)


End[] (* End Private Context *)

EndPackage[]

