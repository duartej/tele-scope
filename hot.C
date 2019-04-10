
// Daniel Pitzl, May 2018
// hot pixels from map

// .x hot.C+("map7")
// .x hot.C+("map0",99)

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
void hot( string hs, int printcut = 999 )
{
  cout << hs;

  TH2I * h2 = (TH2I*)gDirectory->Get( hs.c_str() );

  if( h2 == NULL ) {
    cout << " does not exist in" << gDirectory->GetName() << endl;
    return;
  }

  cout << "  " << h2->GetTitle() << endl;

  // frequency:

  double z9 = 1.1*h2->GetMaximum();
  if( z9 < 1.1 ) z9 = 1.1;

  TH1I * h1 = new
    TH1I( "freq",
	  Form( "%s;counts;pixels",h2->GetZaxis()->GetTitle( ) ),
	  100, 0, z9 );

  const double log10 = log(10);

  TH1I * hl = new
    TH1I( "frel",
	  Form( "%s;log_{10}(counts);pixels",h2->GetZaxis()->GetTitle( ) ),
	  100, 0, log(z9)/log10 );

  //h1->GetXaxis()->SetTitle( h2->GetZaxis()->GetTitle( ) );
  //h1->GetYaxis()->SetTitle( "pixels" );

  multiset <pixel> pxset;

  for( int ii = 1; ii <= h2->GetNbinsX(); ++ii )

    for( int jj = 1; jj <= h2->GetNbinsY(); ++jj ) {

      int n = h2->GetBinContent(ii,jj);

      h1->Fill( n );

      if( n )
	hl->Fill( log( n ) / log10 );
      else
	hl->Fill( -1 ); // underflow = zero

      if( n > 9 ) {
	pixel px{ ii-1, jj-1, n };
	pxset.insert(px);
      }
    }

  cout << "hot " << pxset.size() << endl;

  for( auto px = pxset.begin(); px != pxset.end(); ++px ) {
    cout << "pix "
	 << setw(3) << px->col
	 << setw(5) << px->row
	 << "  " << px->cnt
	 << endl;
    if( px->cnt <= printcut ) break;
  }

  h1->Draw();
  cout << "freq" << endl;

  gPad->SetLogy(1);
  hl->Draw();
  cout << "frel" << endl;

}
