#include <queue>
#include <utility>
#include <vector>
#include <set>
namespace TSnap {

  TFlt GetGroupDegreeCentralization(const int N, TIntH& Degs, const TIntH& GroupNodes, double GD, int k, TIntFltH& mem, TFlt& A, TFlt& B, TIntH& deg_dist, double& approx) {
		TFlt n = N;
		TFlt sum = 0;
		if (k >= 1) {
			k--;
			// fix variables and sum vector
			A *= (n - (double)k) / (n - (double)k - 2.0);
			B *= (n - (double)k - 1) / (n - (double)k - 2.0);
			// iterate sum vector to update values and get the sum
			for (THashKeyDatI<TInt, TFlt> NI = mem.BegI(); NI < mem.EndI(); NI++) {
				int i = NI.GetKey();
				// get old value
				TFlt old_val = NI.GetDat();
				TFlt new_val = old_val * ((n - (double)i - k - 1.0) / (n - (double)k - 2.0));
				// store new value in the sum vector
				mem.AddDat(i, new_val);
				//get the sum
				TFlt deg_val = deg_dist.GetDat(i);
				sum += (deg_val * new_val);
			}
		}
		TFlt centralization = A * GD + sum - B;

		// approximation upper bound
		double approx_rate = 1.58197670686932642438500200510901155854686930107539613626678;
		approx = A * (GD * approx_rate) + sum - B;

		return centralization;
	}

  void DecreaseContribution(TIntIntHH &Histogram, TIntH& Nodes, TIntH& index, int id, int decrease) {
		int HistogramHashId = Nodes.GetDat(id);
		if (Histogram.IsKey(HistogramHashId)) {
			Histogram.GetI(HistogramHashId).GetDat().DelKey(id);
			index.AddDat(HistogramHashId, Histogram.GetI(HistogramHashId).GetDat().Len());
		}
		if (Histogram.IsKey(HistogramHashId - decrease)) {
			Histogram.GetI(HistogramHashId - decrease).GetDat().AddDat(id, id);
			index.AddDat(HistogramHashId - decrease, Histogram.GetI(HistogramHashId - decrease).GetDat().Len());
		}
		else {
			TIntH NewHash;
			NewHash.AddDat(id, id);
			Histogram.AddDat(HistogramHashId - decrease, NewHash);
			index.AddDat(HistogramHashId - decrease, Histogram.GetI(HistogramHashId - decrease).GetDat().Len());
		}
		int prev = Nodes.GetDat(id);
		Nodes.AddDat(id, prev - decrease);
	}

  void MaxCPGreedy(const PNGraph& Graph, TIntH& GroupNodes, TIntFltVH& ResultVector, TFltV& MaxVector, bool& Unimodal, const int k) {
    // buildup container of group nodes
    int br = 0;
    // number of nodes
    int N = Graph->GetNodes();
    // sum of contributions
    double ContributionsSum = 0;
    // Hash table of covered nodes
    TIntH Covered;
    // Indicator if it is in the group set
    TIntBoolH InGroup;
    // DEBUG: file to print out
    //FILE *F = fopen("test10.txt", "wt");
    // id of Result vector element with maximal centralization
    int MaxI = -1;
    // value of maximal centralization
    double MaxGdVal = -1;

    // unimodality check
    bool increasing = true;
    double PrevGdc = -1;
    int NumChanges = 0;

    // Centralization variables
    // sum vector
    TIntFltH mem;
    // degree distribution
    TIntH deg_distribution;
    // variable for updating the sum
    TFlt A = 1 / ((double)N - 1.0);
    // variable for updating the sum
    TFlt B = (double)N / ((double)N - 1.0);

    // Hash table of node contributions (initialy degrees) sorted by values
    TIntH Nodes;
    // Hash table of node in degrees (persuming Degs->GetDat() is faster then NI->GetInDeg())
    TIntH Degs;
    for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      Nodes.AddDat(NI.GetId(), NI.GetInDeg());
      Degs.AddDat(NI.GetId(), NI.GetInDeg());
      InGroup.AddDat(NI.GetId(), false);
    }
    Nodes.SortByDat(false);

    // u of nodes of contributions - Key() - node contribution, Dat() - Hash table of nodes with the distribution 
    TIntIntHH Histogram;
    int PrevDeg = -1;
    TIntH NewHash;
    for (THashKeyDatI<TInt, TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++) {
      int CurrDeg = NI.GetDat();
	  int CurrKey = NI.GetKey();
      if (CurrDeg == PrevDeg) {
        NewHash.AddDat(CurrKey, CurrKey);
      }
      else  {
	    if (PrevDeg != -1) {
          Histogram.AddDat(PrevDeg, NewHash);
		}
        NewHash.Clr();
        NewHash.AddDat(CurrKey, CurrKey);
        PrevDeg = CurrDeg;
      }
    }
    if (NewHash.Len() > 0) {
      Histogram.AddDat(PrevDeg, NewHash);
    }
    Histogram.SortByKey(false);

    // Index of histogram values - index is acctualy a histogram of contributions. Key() - contribution, Dat() - number of nodes with that contribution
    TIntH index;
    for (THashKeyDatI<TInt, TIntH> HI = Histogram.BegI(); HI < Histogram.EndI(); HI++) {
      index.AddDat(HI.GetKey(), HI.GetDat().Len());
      mem.AddDat(HI.GetKey(), 1 / ((double)N - 1));
      deg_distribution.AddDat(HI.GetKey(), HI.GetDat().Len());
    }

    // Choose node
    int ChoosenNode = -1;
    // Choose the biggest degree
    int ChoosenDeg = Histogram.BegI().GetKey();
    // Group nodes counter
    int GroupNodesCounter = 0;

    // Main loop
    do {
      // DEBUG: Print out the graph
      // printf("\nGraph:"); for (TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) { printf("\n %u -> %i", EI.GetSrcNId(), EI.GetDstNId()); }
      // DEBUG: Print out the contributions
      // printf("\nContributions:"); for (THashKeyDatI<TInt, TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++) { printf("\n %i: %i", NI.GetKey(), NI.GetDat()); }
      // DEBUG Histogram:
      // printf("\nHistogram:"); for (THashKeyDatI<TInt, TIntH> HI = Histogram.BegI(); HI < Histogram.EndI(); HI++) { TIntH deg = HI.GetDat(); printf("\n %i(%i):", HI.GetKey(), deg.Len()); for (THashKeyDatI<TInt, TInt> NI = deg.BegI(); NI < deg.EndI(); NI++){ printf("%i,", NI.GetDat()); } }
      // THashKeyDatI<TInt, TInt> FirstIndex = index.BegI();
      // Reset the choosen node
      // Main loop - start timer
      double tm1, tm2;
      TExeTm ExeTm;
      tm1 = ExeTm.GetSecs();
      ChoosenNode = -1;
      // Choosing the node - this takes in sum O(maxdeg) steps, or O(n) steps
      while (ChoosenNode == -1 && ChoosenDeg >= 0) {
        if (index.IsKey(ChoosenDeg)) {
          if (index.GetDat(ChoosenDeg) > 0)
            ChoosenNode = ChoosenNode = Histogram.GetDat(ChoosenDeg).BegI().GetDat();
          else
            ChoosenDeg--;
        }
        else {
          ChoosenDeg--;
        }
      }

      // Choosing node end timer
      double tm4 = ExeTm.GetSecs();
      double time4 = tm4 - tm1;
      if (ChoosenNode != -1) {
        ContributionsSum += ChoosenDeg;
        GroupNodes.AddDat(GroupNodesCounter, ChoosenNode);
        GroupNodesCounter++;
        // If the choosen node is in the set of covered neighbours it has to be removed
        if (Covered.IsKey(ChoosenNode)) {
          Covered.DelKey(ChoosenNode);
        }
        // get all the nodes with choosen degree
        TIntH FirstHash = Histogram.GetDat(ChoosenDeg);
        // first key = choosen deg
        int FirstKey = ChoosenDeg;
        // Number of nodes with choosen degree
        int Len = FirstHash.Len();

        // Remove the index for this degree if it is the last element, otherwise just decrease the size
        if (Len == 0 || Len == 1) {
          Histogram.DelKey(ChoosenDeg);
          index.DelKey(ChoosenDeg);
        }
        else {
          Histogram.GetI(ChoosenDeg).GetDat().DelKey(ChoosenNode);
          index.AddDat(FirstKey, Len - 1);
        }

        // Get the node
        TNGraph::TNodeI NI1 = Graph->GetNI(ChoosenNode);

        // DEBUG:
        int outdeg = NI1.GetOutDeg();

        // A) reduce contributionto in-degree nodes by 1 (or 2)
        for (int i = 0; i < NI1.GetInDeg(); i++) {
          int nei = NI1.GetInNId(i);
          // DEBUG:
          // printf("\n %i <- %i", ChoosenNode, nei);
          int decrease = 1;
          if (Nodes.IsKey(nei) && ChoosenNode != nei) {
            DecreaseContribution(Histogram, Nodes, index, nei, decrease);
          }
        }

        // B) for all node-outneigh neigbours
        for (int i = 0; i < NI1.GetOutDeg(); i++) {
	      // Id of out neighbout
          int n_out = NI1.GetOutNId(i);
          DecreaseContribution(Histogram, Nodes, index, n_out, 1);
	      if (!Covered.IsKey(n_out)) {
		     Covered.AddKey(n_out);
	      }
          if (ChoosenNode != n_out) {
          // DEBUG:
		  // printf("\n %i -> %i", ChoosenNode, n_out);

          // Get all node-outneigh-inneigh nodes
          TNGraph::TNodeI NI1 = Graph->GetNI(n_out);
          TIntV n_out_ins;
          for (int j = 0; j < NI1.GetInDeg(); j++) {
            // Id of node-outneigh-inneigh
            int n_out_in = NI1.GetInNId(j);
            if (n_out_in != ChoosenNode) {
              n_out_ins.Add(n_out_in);
            }
          }

          if (n_out_ins.Len() > 0) {
            for (int i = 0; i < n_out_ins.Len(); i++) {
              int n_out_in = n_out_ins[i];
              // DEBUG:
              // printf("\n %i -> %i <- %i", ChoosenNode, n_out, n_out_in);
              int HistogramHashId = Nodes.GetDat(n_out_in);

              // B.1 Delete in edges
              if (Graph->IsEdge(n_out_in, n_out)) {
                Graph->DelEdge(n_out_in, n_out);
              }

              // B.2 DECREASE node-outneigh-inneigh for one
              DecreaseContribution(Histogram, Nodes, index, n_out_in, 1);
              // DEBUG:
              // printf("\n-%i", n_out_in);
            }
          }
        }
      }

      // C) delete the node
      Graph->DelNode(ChoosenNode);
      Nodes.DelKey(ChoosenNode);

      // Main Loop end timer
      double tm2 = ExeTm.GetSecs();
      double time = tm2 - tm1;

      // Timer for CENTRALIZATION
      // DEBUG:
      // for (int i = 0; i < Covered.Len(); i++) { printf("\n %i. %i", i + 1, Covered[i]); }
      // DEBUG: ContributionsSum should be == to Covered.Len()
      double approx = 0;
      TFlt gdc2 = TSnap::GetGroupDegreeCentralization(N, Degs, GroupNodes, /*Covered.Len()*/ContributionsSum, GroupNodesCounter, mem, A, B, deg_distribution, approx);
      double tm3 = ExeTm.GetSecs();
      double time1 = tm3 - tm2;

      // PRINT OUT
      // printf("\n %i %i %i %i %f %f %i %f %f %f %f", GroupNodesCounter, ChoosenNode, ChoosenDeg, Degs.GetDat(ChoosenNode), ContributionsSum, ContributionsSum + GroupNodesCounter, N, gdc2, time, time1, time4);
      if (GroupNodesCounter % 1000 == 0) {
        printf("\n %f %", (ContributionsSum + (double)GroupNodesCounter) / (double)N);
      }
      // DEBUG:
      //fprintf(F, "\n %i %i %i %i %f %f %i %f %f %f %f", GroupNodesCounter, ChoosenNode, ChoosenDeg, Degs.GetDat(ChoosenNode), ContributionsSum, ContributionsSum + GroupNodesCounter, N, gdc2, time, time1, time4);

      // save result vector
      TFltV vec;
      vec.Add((double)GroupNodesCounter); vec.Add(ChoosenNode); vec.Add(ChoosenDeg); vec.Add((double)Degs.GetDat(ChoosenNode)); vec.Add(ContributionsSum); vec.Add(gdc2); vec.Add(approx);
      ResultVector.AddDat(GroupNodesCounter, vec);

      // check if its max
      if (gdc2 > MaxGdVal) {
        MaxI = GroupNodesCounter;
        MaxGdVal = gdc2;
        MaxVector = vec;
      }

        // count number of changes of unimodality trend
        if (PrevGdc <= gdc2 && !increasing) {
          NumChanges++;
          increasing = true;
        }
        else if (PrevGdc >= gdc2 && increasing) {
          NumChanges++;
          increasing = false;
        }
        else {}
        PrevGdc = gdc2;
      }

    } while (ChoosenNode > -1);

    if (NumChanges > 1) {
      Unimodal = false;
    }
  }

}; // namespace TSnap