(* ::Package:: *)

(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Dec 27, 2021 *)

BeginPackage["chiralGraph`"]
(* Exported symbols added here with SymbolName::usage *) 

negativeAndZMspectrum::usage = "negativeAndZMspectrum[hamiltonian] returns {{negEval,negEvec},{zmEval,zmEvec}}, where negEval and zmEval are the negative and zero eigenvalues, and negEvec and zmEvec are the negative and zero eigenvectors, respectively."
negativeAndZMspectrum::badarg = "The function expects a bidimensional list as an argument."

localizedBasisPencilMatrix::usage = "localizedBasisPencilMatrix[graph,basis] gives {localisedBasis, accuracyList}: the basis of localised function obtained with the pencil matrix method for the given graph and given negative eigenvalues and eigenvectors, and the accuracy of the eigenvectors."
localizedBasisPencilMatrix::badarg = "The function expects a chiral graph, a list og negative eigenvalues, and a list of the associated eigenvectors."

localizedBasisPencilMatrix3D::usage = "localizedBasisPencilMatrix[graph,basis] gives {localisedBasis, accuracyList}: the basis of localised function obtained with the pencil matrix method for the given graph and given negative eigenvalues and eigenvectors, and the accuracy of the eigenvectors."
localizedBasisPencilMatrix3D::badarg = "The function expects a chiral graph, a list og negative eigenvalues, and a list of the associated eigenvectors."

localizedWannier::usage = "localizedBasisPencilMatrix[graph,basis,nIterations:200,step:0.1] gives the basis of maximally localized Wannier functions for the given graph and given basis."
localizedWannier::badarg = "The function expects a chiral graph and a basis."

localizedWannier3D::usage = "localizedBasisPencilMatrix[graph,basis,nIterations:200,step:0.1] gives the basis of maximally localized Wannier functions for the given graph and given basis."
localizedWannier3D::badarg = "The function expects a chiral graph and a basis."

chiralPolarizationField::usage = "chiralPolarizationField[graph,localBasis] gives the chiral polarization field {x_list,y_list,Pix_list,Piy_list} the graph using localizedBasis."
chiralPolarizationField::badarg = "The function expects a graph and a list containing the localized basis as arguments."

chiralPolarizationField3D::usage = "chiralPolarizationField[graph,localBasis] gives the chiral polarization field {x_list,y_list,z_list,Pix_list,Piy_list,Piz_list} the graph using localizedBasis."
chiralPolarizationField3D::badarg = "The function expects a graph and a list containing the localized basis as arguments."

plotChiralPolarizaionField::usage = "plotChiralPolarizaionField[{x_list,y_list,Pix_list,Piy_list}] returns a graphic object with the chiral polarization field, and a graphic object with the color disk."
plotChiralPolarizaionField::badarg = "The function expects a list of the format {x_list,y_list,Pix_list,Piy_list}."

chiralList::usage = "chiralList[chiralGraph] gives a list ordered by VertexList assigning +1 to dofs and -1 to restrictions"
chiralList::badarg = "wrong format/arguments"

projAList::usage = "projAList[chiralGraph] gives a list ordered by VertexList assigning +1 to dofs and 0 to restrictions"
projAList::badarg = "wrong format/arguments"

projBList::usage = "projBList[chiralGraph] gives a list ordered by VertexList assigning 0 to dofs and +1 to restrictions"
projBList::badarg = "wrong format/arguments"

positionList::usage ="prositionList[chiralGraph] gives a list ordered by Vertexlist assigning the respective positions."
positionList::badarg = "wrong format/arguments"

highlightGraph::usage ="highlightGraph[graph,state,opacity] highlights the given state on the given graph"
highlightGraph::badarg = "wrong format/arguments";

unitaryEvolution::usage = "unitaryEvolution[hamiltonian,initialState,{ti,tf,dt}] returns a list with all the states from ti to tf with time steps dt and unitary operator from the hamiltonian."
unitaryEvolution::badarg = "wrong format/arguments";

spreading::usage = "spreading[positions_List,basis_List] returns {spread,\[Del]spread} of the given basis"
spreading::badarg = "wrong format/arguments"

Begin["`Private`"]
(* Implementation of the package *)

negativeAndZMspectrum[hamiltonian_List/;Length@Dimensions[hamiltonian]==2,threshold:(_?Positive):0.03,offsetEnergy_:0] :=
    Module[ {sortEigenVal,sortEigenVec},
        {sortEigenVal, sortEigenVec} = 
        SortBy[Transpose[N@Eigensystem[hamiltonian]], First] // Transpose;
        Module[ {f1, f2, eigenSelect},
            f1 = # < offsetEnergy - threshold &;
            f2 = Abs[offsetEnergy-#] < threshold &;
            eigenSelect[eigenval_, eigenvec_, 
              booleanFunction_] :=
                {Select[eigenval, booleanFunction], 
                Pick[eigenvec, booleanFunction[#] & /@ eigenval]};
            eigenSelect[sortEigenVal, sortEigenVec, #] & /@ {f1, f2}
        ]
    ]
negativeAndZMspectrum[___] :=
    (Message[negativeAndZMspectrum::badarg];
     $Failed)

chiralList[graph_Graph] :=
    With[ {rgbToColor = If[ #===RGBColor[0, 0, 1],
                            -1,
                            1
                        ]&},
        rgbToColor[AnnotationValue[{graph, #}, VertexStyle]] & /@ 
        VertexList[graph]
    ]
chiralList[___] :=
    (Message[chiralList::badarg];
     $Failed)

projAList[graph_Graph] :=
	With[ {rgbToColor = If[ #===RGBColor[0, 0, 1],
                            0,
                            1
                        ]&},
        rgbToColor[AnnotationValue[{graph, #}, VertexStyle]] & /@ 
        VertexList[graph]
    ];
projAList[___] :=
    (Message[projAList::badarg];
     $Failed)    
     
projBList[graph_Graph] :=
	With[ {rgbToColor = If[ #===RGBColor[0, 0, 1],
                            1,
                            0
                        ]&},
        rgbToColor[AnnotationValue[{graph, #}, VertexStyle]] & /@ 
        VertexList[graph]
    ];
projBList[___] :=
    (Message[projBList::badarg];
     $Failed)  


positionList[graph_Graph] :=
    AnnotationValue[graph, VertexCoordinates]
positionList[___] :=
    (Message[positionList::badarg];
     $Failed)  

matrixPencil[m1_List,m2_List,l_]:= l m1+(1-l) m2
matrixPencil[m1_List,m2_List,{l1_,l2_}]:= l1 m1+l2 m2
matrixPencil[m1_List,m2_List,m3_List,{l1_,l2_}]:= l1 m1+l2 m2+(1-l1-l2) m3
matrixPencil[___] := $Failed


localizedBasisPencilMatrix[graph_Graph,negEvec_List,l_:0.3] :=
    Module[ {localFunctions,projections, projx, projy, sortPMEval,
    sortPMEvec,positions},
        projections[posComponent_] :=
            Outer[Conjugate[#1] . (posComponent #2) &, negEvec, negEvec, 1];
        positions=positionList[graph];
        {projx, projy} = projections[#] & /@ (positions\[Transpose]);
        {sortPMEval, sortPMEvec} = 
         SortBy[Eigensystem[matrixPencil[projx,projy,l]]\[Transpose],
            First]\[Transpose];
        localFunctions = sortPMEvec . negEvec;
        {localFunctions, spreading[positions, localFunctions][[1]]}
    ]
localizedBasisPencilMatrix[graph_Graph,negEvec_List,{l1_:0.3,l2_:0.7}] :=
    Module[ {localFunctions,projections, projx, projy, sortPMEval,
    sortPMEvec,positions},
        projections[posComponent_] :=
            Outer[Conjugate[#1] . (posComponent #2) &, negEvec, negEvec, 1];
        positions=positionList[graph];
        {projx, projy} = projections[#] & /@ (positions\[Transpose]);
        {sortPMEval, sortPMEvec} = 
         SortBy[Eigensystem[matrixPencil[projx,projy,{l1,l2}]]\[Transpose],
            First]\[Transpose];
        localFunctions = sortPMEvec . negEvec;
        {localFunctions, spreading[positions, localFunctions][[1]]}
    ]
localizedBasisPencilMatrix[___] := (Message[localizedBasisPencilMatrix::badarg];$Failed)


localizedBasisPencilMatrix3D[graph_Graph,negEvec_List,{l1_:0.3,l2_:0.7}] :=
    Module[ {localFunctions,projections, projx, projy, projz, sortPMEval,
    sortPMEvec,positions},
        projections[posComponent_] :=
            Outer[Conjugate[#1] . (posComponent #2) &, negEvec, negEvec, 1];
        positions=positionList[graph];
        {projx, projy, projz} = projections[#] & /@ (positions\[Transpose]);
        {sortPMEval, sortPMEvec} = 
         SortBy[Eigensystem[matrixPencil[projx,projy,projz,{l1,l2}]]\[Transpose],
            First]\[Transpose];
        localFunctions = sortPMEvec . negEvec;
        {localFunctions, spreading3D[positions, localFunctions][[1]]}
    ]
localizedBasisPencilMatrix3D[___] := (Message[localizedBasisPencilMatrix3D::badarg];$Failed)


spreading[positions_List, basis_List] := 
 Module[{projCX, projCY, xD, yD, xP, yP, nNeg,gradient},
  nNeg = Length@basis; 
  projCX = 
   Chop[Table[
    Conjugate[basis[[n]]] . (positions[[All, 1]] basis[[m]]), {n, 1, nNeg}, {m, 
     1, nNeg}],10^-7];
  projCY = 
   Chop[Table[
    Conjugate[basis[[n]]] . (positions[[All, 2]] basis[[m]]), {n, 1, nNeg}, {m, 
     1, nNeg}],10^-7];
  {xD, yD} = DiagonalMatrix[Diagonal[#]] & /@ {projCX, projCY};
  {xP, yP} = {projCX - xD, projCY - yD};
  gradient = Chop[2 (xP . xD - xD . xP + yP . yD - yD . yP),10^-7];
  (*Nor normalized by the size of UC*)
  {Chop[Tr[xP . xP + yP . yP],10^-7],gradient}
  ]
spreading[___] := (Message[spreading::badarg];$Failed)


spreading3D[positions_List, basis_List] := 
 Module[{projCX, projCY, projCZ, xD, yD, zD, xP, yP, zP, nNeg,gradient},
  nNeg = Length@basis; 
  projCX = 
   Chop@Table[
    basis[[n]] . (positions[[All, 1]] basis[[m]]), {n, 1, nNeg}, {m, 
     1, nNeg}];
  projCY = 
   Chop@Table[
    basis[[n]] . (positions[[All, 2]] basis[[m]]), {n, 1, nNeg}, {m, 
     1, nNeg}];
  projCZ = 
   Chop@Table[
    basis[[n]] . (positions[[All, 3]] basis[[m]]), {n, 1, nNeg}, {m, 
     1, nNeg}];
  {xD, yD, zD} = DiagonalMatrix[Diagonal[#]] & /@ {projCX, projCY, projCZ};
  {xP, yP, zP} = {projCX - xD, projCY - yD, projCZ-zD};
  gradient = 2 (xP . xD - xD . xP + yP . yD - yD . yP + zP . zD - zD . zP);
  (*Nor normalized by the size of UC*)
  {Tr[xP . xP + yP . yP + zP . zP],gradient}
  ]
spreading3D[___] := (Message[spreading3D::badarg];$Failed)



localizedWannier[graph_Graph,basis_List,nIterations_:200,stepSize_:0.001] :=
    Module[ {functionalTable,xD,yD,xP,yP,projCX,projCY,positions,gradient,localFunctions,nNeg,spread},
        localFunctions = basis;
        nNeg = Length@basis;
        positions = positionList[graph];
        functionalTable = 
          Table[
          	(*projCX = 
            Table[
             localFunctions[[n]] . (positions[[All, 1]] localFunctions[[m]]), {n, 1, 
              nNeg}, {m, 1, nNeg}];
                
            projCY = 
                 Table[
                  localFunctions[[n]] . (positions[[All, 2]] localFunctions[[m]]), {n, 1, 
                   nNeg}, {m, 1, nNeg}];
                {xD, yD} = DiagonalMatrix[Diagonal[#]] & /@ {projCX, projCY};
                {xP, yP} = {projCX - xD, projCY - yD};
                gradient = 2 (xP . xD - xD . xP + yP . yD - yD . yP);*)
                {spread,gradient} = spreading[positions, localFunctions];
                localFunctions = #/Norm[#] & /@ (localFunctions - stepSize gradient . localFunctions);
                (*spread not normalized by the size of the unit cell*)
                spread, nIterations
           ];
           {localFunctions,functionalTable}
    ]
localizedWannier[___] := (Message[localizedWannier::badarg];$Failed)


localizedWannier3D[graph_Graph,basis_List,nIterations_:200,stepSize_:0.001] :=
    Module[ {functionalTable,xD,yD,xP,yP,projCX,projCY,positions,gradient,localFunctions,nNeg,spread},
        localFunctions = basis;
        nNeg = Length@basis;
        positions = positionList[graph];
        functionalTable = 
          Table[
          	(*projCX = 
            Table[
             localFunctions[[n]] . (positions[[All, 1]] localFunctions[[m]]), {n, 1, 
              nNeg}, {m, 1, nNeg}];
                
            projCY = 
                 Table[
                  localFunctions[[n]] . (positions[[All, 2]] localFunctions[[m]]), {n, 1, 
                   nNeg}, {m, 1, nNeg}];
                {xD, yD} = DiagonalMatrix[Diagonal[#]] & /@ {projCX, projCY};
                {xP, yP} = {projCX - xD, projCY - yD};
                gradient = 2 (xP . xD - xD . xP + yP . yD - yD . yP);*)
                {spread,gradient} = spreading3D[positions, localFunctions];
                localFunctions = #/Norm[#] & /@ (localFunctions - stepSize gradient . localFunctions);
                (*spread not normalized by the size of the unit cell*)
                spread, nIterations
           ];
           {localFunctions,functionalTable}
    ]
localizedWannier3D[___] := (Message[localizedWannier::badarg];$Failed)

chiralPolarizationField[graph_Graph,localBasis_List] := With[{positions=positionList[graph],chiral=chiralList[graph]},Transpose[({Conjugate[#] . (positions[[All, 
          1]] #), Conjugate[#] . (positions[[All, 2]] #), 
      2 Conjugate[#] . (positions[[All, 1]] chiral #), 
      2 Conjugate[#] . (positions[[All, 2]] chiral #)} & /@ 
    localBasis)]]   
chiralPolarizationField[___] := (Message[chiralPolarizationField::badarg];$Failed) 


chiralPolarizationField3D[graph_Graph,localBasis_List] := With[{positions=positionList[graph],chiral=chiralList[graph]},
	Transpose[({
		Conjugate[#] . (positions[[All, 1]] #),
		Conjugate[#] . (positions[[All, 2]] #),
		Conjugate[#] . (positions[[All, 3]] #),
		2 Conjugate[#] . (positions[[All, 1]] chiral #),
		2 Conjugate[#] . (positions[[All, 2]] chiral #),
		2 Conjugate[#] . (positions[[All, 3]] chiral #)} & /@ 
    localBasis)]]   
chiralPolarizationField3D[___] := (Message[chiralPolarizationField3D::badarg];$Failed)  
     
plotChiralPolarizaionField[{x_List,y_List,piX_List,piY_List},fracArrowHead_:0.038,absThickness_:2] :=
    Module[ {customColorF,colorDisk,colorchiralPolarizationField},
    	customColorF = Blend[{RGBColor[0.685695, 0.242449, 0.268261], RGBColor[
    0.534081, 0.0853132, 0.16669], GrayLevel[0], RGBColor[
    0.139681, 0.311666, 0.550652], RGBColor[
    0.256859, 0.523007, 0.711644], RGBColor[
    0.433786, 0.670834, 0.793785], RGBColor[0.289326, 0.685107, 0.24], 
    RGBColor[0., 0.442859, 0.0749256], RGBColor[
    0.44, 0.685107, 0.108759], RGBColor[0.87, 0.67, 0.522424], 
    RGBColor[0.685695, 0.242449, 0.268261]}, #] &;
        colorDisk = 
        With[ {offset = -1.5, sectors = 360., angleAux = 2. Pi/360},
            Graphics[
             Table[{customColorF[i/sectors], 
               EdgeForm[{Thick, customColorF[i/sectors]}], 
               Annulus[{0, 0}, {1/2, 
                 1}, {i angleAux + offset, (i + 1) angleAux + offset}]}, {i, 0,
                sectors - 1}], ImageSize -> Tiny]
        ];
        colorchiralPolarizationField=With[ {enlargementFactor = 1, 
        angles = MapThread[ArcTan[#1, #2] &, {piX, piY}]},
            Graphics[
             MapThread[{customColorF[Mod[1 + #5/(2 \[Pi]) + 0.2, 1]], 
                Arrowheads[fracArrowHead Norm[{#3, #4}]], AbsoluteThickness[absThickness], 
                Arrow[{{#1 - enlargementFactor #3/2, #2 - 
                    enlargementFactor #4/2}, {#1 + enlargementFactor #3/2, #2 + 
                    enlargementFactor #4/2}}]} &, {x, 
                y, piX, piY}~Join~{angles}]]
        ];
        {colorchiralPolarizationField,colorDisk}
    ]
plotChiralPolarizaionField[{x_List,y_List,z_List,piX_List,piY_List,piZ_List},fracArrowHead_:0.038,absThickness_:2] :=
    Module[ {customColorF,colorDisk,colorchiralPolarizationField},
    	customColorF = Blend[{RGBColor[0.685695, 0.242449, 0.268261], RGBColor[
    0.534081, 0.0853132, 0.16669], GrayLevel[0], RGBColor[
    0.139681, 0.311666, 0.550652], RGBColor[
    0.256859, 0.523007, 0.711644], RGBColor[
    0.433786, 0.670834, 0.793785], RGBColor[0.289326, 0.685107, 0.24], 
    RGBColor[0., 0.442859, 0.0749256], RGBColor[
    0.44, 0.685107, 0.108759], RGBColor[0.87, 0.67, 0.522424], 
    RGBColor[0.685695, 0.242449, 0.268261]}, #] &;
        colorDisk = 
        With[ {offset = -1.5, sectors = 360., angleAux = 2. Pi/360},
            Graphics[
             Table[{customColorF[i/sectors], 
               EdgeForm[{Thick, customColorF[i/sectors]}], 
               Annulus[{0, 0}, {1/2, 
                 1}, {i angleAux + offset, (i + 1) angleAux + offset}]}, {i, 0,
                sectors - 1}], ImageSize -> Tiny]
        ];
        colorchiralPolarizationField=With[ {enlargementFactor = 1},
            Graphics3D[
             MapThread[{Arrowheads[fracArrowHead Norm[{#4, #5, #6}]], AbsoluteThickness[absThickness], 
                Arrow[{{#1 - enlargementFactor #4/2, #2 - 
                    enlargementFactor #5/2, #3 - 
                    enlargementFactor #6/2}, {#1 + enlargementFactor #4/2, #2 + 
                    enlargementFactor #5/2,  #3 + 
                    enlargementFactor #6/2}}]} &, {x, 
                y, z, piX, piY, piZ}]]
        ];
        {colorchiralPolarizationField,colorDisk}
    ]
plotChiralPolarizaionField[___] :=
    (Message[plotChiralPolarizaionField::badarg];
     $Failed)


(*Highlight function*)
Clear[highlightGraph];
highlightGraph[graph_Graph,state_List,baseOpacity_:0.2]:=With[{cList=chiralList[graph]},Graph[graph,VertexStyle->MapThread[#1->Directive[EdgeForm[None],If[#3==1,Red,Blue],Opacity[baseOpacity+#2]]&,{VertexList[graph],state,cList}]]];
highlightGraph[___]:=(Message[highlightGraph::badarg];$Failed)


unitaryEvolution[hamiltonian_, 
   initialState_, {ti_, 
    tf_, \[CapitalDelta]t_}] /; (Length[initialState] == 
    Length[hamiltonian]) := Module[{tRange},
  tRange = Range[ti, tf, \[CapitalDelta]t];
  With[{uu = N@MatrixExp[-I hamiltonian \[CapitalDelta]t]}, 
   NestList[uu . # &, initialState, Length@tRange - 1]]]
unitaryEvolution[___]:=(Message[unitaryEvolution::badarg];
     $Failed)

End[]

EndPackage[]

