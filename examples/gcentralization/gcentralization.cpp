#include "stdafx.h"
int main(int argc, char* argv[]) {
  const double M_e = 2.7182818284590452353602874;
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Group Centralization. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm, ExeTm1, ExeTm2;
  double tm1 = 0, tm2 = 0;
  Try

	// loading input parameters
    const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "./Data/cobiss.txt", "Input un/directed graph");
    const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "./Output/gc.out.txt", "Output file");

    // group nodes for output
    TIntH gk1;
    // result vectors for oputput
    TIntFltVH ResultVector;
    // max vector for output
    TFltV MaxVector;
    // unimodality variable for output
    bool Unimodal = true;

    // loading graph
    printf("Loading %s...\n", InFNm.CStr());
    PUNGraph UGraph = TSnap::LoadEdgeList<PUNGraph>(InFNm);
    // printing the size of node and edge lists
    printf("nodes:%d  edges:%d\n", UGraph->GetNodes(), UGraph->GetEdges());

    // Preprocesing - creating a proper directed graph without loops
    printf("Preprocessing %s...\n", InFNm.CStr());
    // converting to directed graph
    //PNGraph DGraph = TSnap::ConvertGraph<PNGraph, PUNGraph>(UGraph);
    PNGraph DGraph = TNGraph::New();
    // generating directed graph

    for (TUNGraph::TNodeI NI = UGraph->BegNI(); NI < UGraph->EndNI(); NI++) {
      DGraph->AddNode(NI.GetId());
    }

    for (TUNGraph::TEdgeI EI = UGraph->BegEI(); EI < UGraph->EndEI(); EI++) {
      int n1 = EI.GetSrcNId();
      int n2 = EI.GetDstNId();
      if (n1 != n2) { // leaving out loops
        DGraph->AddEdge(n1, n2);
        DGraph->AddEdge(n2, n1);
      }
	}

    // remove loops - if we dont create a new directed graph, but use the initial graph, we have to make sure the loops are removed
	/*
    for (TNGraph::TNodeI NI = DGraph->BegNI(); NI < DGraph->EndNI(); NI++) {
      if (DGraph->IsEdge(NI.GetId(), NI.GetId())) {
        DGraph->DelEdge(NI.GetId(), NI.GetId());
      }
    }
	*/

    // end preprocessing

    // max size of the group is initialy the number of nodes
    int k = UGraph->GetNodes();

    // Main part - centralization
    printf("Centralization...\n");
    tm1 = ExeTm1.GetSecs();
    TSnap::MaxCPGreedy(DGraph, gk1, ResultVector, MaxVector, Unimodal, k);
    tm2 = ExeTm2.GetSecs();
    printf("\nExecution time %f\n", tm2 - tm1);

    // generating plot
    TFltPrV prv, prv1, prv2;

    // print results to the file
    printf("Printing results...\n");
    FILE *F = fopen(OutFNm.CStr(), "wt");
	for (THashKeyDatI<TInt, TFltV> VI = ResultVector.BegI(); VI < ResultVector.EndI(); VI++) {
      TFltV vec = VI.GetDat();
      int id = VI.GetKey();
      fprintf(F, "%i ", id);
      for (int i = 0; i < vec.Len(); i++) {
        fprintf(F, "%f ", vec[i]);
      }

    // centralization plot data
    TFltPr p;
    p.Val1 = vec[0]; p.Val2 = vec[5]; prv.Add(p);

    // centralization app border plot data
    TFltPr p2;
    p2.Val1 = vec[0]; p2.Val2 = vec[6]; prv2.Add(p2);

    // centrality plot data
    TFltPr p1;
    p1.Val1 = vec[0]; p1.Val2 = vec[4]; prv1.Add(p1);

    fprintf(F, "\n");
  }

  // print the max result
  fprintf(F, "***\n");
  for (int i = 0; i < MaxVector.Len(); i++) {
    fprintf(F, "%f ", MaxVector[i]);
  }
  fprintf(F, "\nUnimodal: %s", Unimodal ? "true" : "false");
  fprintf(F, "\n***");

  // printing plots
  TGnuPlot::PlotValV(prv, TStr::Fmt("%s-centralization-k", OutFNm.CStr()), "Centralization to group size", "Group size", "Centralization", gpsAuto, false, gpwLinesPoints);
  TGnuPlot::PlotValV(prv1, TStr::Fmt("%s-centrality-k", OutFNm.CStr()), "Centrality to group size", "Group size", "Centralization", gpsAuto, false, gpwLinesPoints);
  
  // plot with the bound
  /*TGnuPlot GnuPlot(TStr::Fmt("%s-centralization1-k", OutFNm.CStr()), "Centralization to group size", false);
  GnuPlot.AddPlot(prv, gpwLinesPoints, "Centralization", "pt 8");
  GnuPlot.AddPlot(prv2, gpwLines, "Upper bound", "pt 8");
  GnuPlot.SetXYLabel("Group size", "Centralization");
  GnuPlot.SetScale(gpsAuto);
  GnuPlot.SavePng();*/

  printf("Done!\n");
  system("pause");

  Catch
    return 0;
}