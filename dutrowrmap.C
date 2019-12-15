
// D. Pitzl, Nov 2018
// row correlation coefficient map for scope53

// .x dutrowrmap.C+

#include <TDirectory.h>
#include <TProfile2D.h>
#include <iostream>
using namespace std;

void dutrowrmap()
{
  TProfile2D * pq = (TProfile2D*)gDirectory->Get( "dutrowpqmap" );
  TProfile2D * pp = (TProfile2D*)gDirectory->Get( "dutrowppmap" );
  TProfile2D * qq = (TProfile2D*)gDirectory->Get( "dutrowqqmap" );

  if( pq == NULL ) {
    cout << "pq does not exist\n";
    return;
  }
  if( pp == NULL ) {
    cout << "pp does not exist\n";
    return;
  }
  if( qq == NULL ) {
    cout << "qq does not exist\n";
    return;
  }

  TProfile2D * rr = (TProfile2D*)pq->Clone();
  rr->Reset();
  rr->SetTitle("DUT previous row correlation");
  rr->GetZaxis()->SetTitle("DUT previous row signal correlation");

  for( int i = 1; i <= pq->GetNbinsX(); ++i ) {
    cout << " " << i << flush;
    for( int j = 1; j <= pq->GetNbinsY(); ++j ) {
      int ij = pq->GetBin(i,j);
      if( pq->GetBinEntries(ij) > 2 && pp->GetBinContent(i,j) * qq->GetBinContent(i,j) > 0 ) {
	double r = pq->GetBinContent(i,j) / sqrt( pp->GetBinContent(i,j) * qq->GetBinContent(i,j) );
	rr->Fill( pq->GetXaxis()->GetBinCenter(i), pq->GetXaxis()->GetBinCenter(j), r );
	/*
	cout << i << " " << j  << " " << ij 
	     << " " << pq->GetBinEntries(ij)
	     << " " << pq->GetBinContent(i,j)
	     << " " << pp->GetBinContent(i,j)
	     << " " << qq->GetBinContent(i,j)
	     << " " << r
	     << endl;
	*/
      }
    }
  }
  cout << endl;
  rr->Draw("colz");
}
