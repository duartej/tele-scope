
// Daniel Pitzl, May 2018
// hot pixels from map

// .x hot.C+("map7")
// .x hot.C+("map8",9999)  RD53

#include "TDirectory.h"
#include "TH2I.h"
#include "TPad.h"
#include <set> // multiset is in set
#include <iostream> // cout
#include <iomanip> // setw

struct pixel
{
  int col;
  int row;
  int cnt;
  bool operator < (const pixel & pxObj ) const
  {
    return cnt > pxObj.cnt;
  }
};

//------------------------------------------------------------------------------
void hot( string hs, int hotcut = 999 )
{
  cout << hs;

  TH2I * h2 = (TH2I*)gDirectory->Get( hs.c_str() );

  if( h2 == NULL ) {
    cout << " does not exist in" << gDirectory->GetName() << endl;
    return;
  }

  cout << "  " << h2->GetTitle() << endl;

  // frequency:

  double z9 = 1.02*h2->GetMaximum();
  if( z9 < 1.1 ) z9 = 1.1;

  TH1I * h1 = new
    TH1I( "freq",
	  Form( "%s;hits/pixel;pixels",h2->GetZaxis()->GetTitle( ) ),
	  200, 0, z9 );

  const double log10 = log(10);

  TH1I * hl = new
    TH1I( "frel",
	  Form( "%s;log_{10}(hits/pixel);pixels",h2->GetZaxis()->GetTitle( ) ),
	  100, 0, log(z9)/log10 );

  multiset <pixel> pxset;

  for( int ii = 1; ii <= h2->GetNbinsX(); ++ii )

    for( int jj = 1; jj <= h2->GetNbinsY(); ++jj ) {

      int n = ( h2->GetBinContent(ii,jj) + 0.1 );

      if( n ) { // active
	h1->Fill( n );
	hl->Fill( log( n ) / log10 );

	pixel px{ ii-1, jj-1, n };
	pxset.insert(px);
      }

    }

  int nhot = 0;
  for( auto px = pxset.begin(); px != pxset.end(); ++px ) {
    cout << "pix "
	 << setw(3) << px->col
	 << setw(5) << px->row
	 << "  " << px->cnt
	 << endl;
    ++nhot;
    if( px->cnt <= hotcut ) break;
  }

  cout << "hot " << nhot << " above " << hotcut << endl;

  cout << "active " << pxset.size() << endl;

  h1->Draw();
  cout << "freq" << endl;

  gPad->SetLogy(1);
  hl->Draw();
  cout << "frel" << endl;

}
