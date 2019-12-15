
// Daniel Pitzl, Nov 2019
// RFD53A thr scan
// root -l thr563.C

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>

//------------------------------------------------------------------------------
void thr563()
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleOffset( 1.4, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(111111);
  gStyle->SetStatFormat("8.6g"); // more digits, default is 6.4g
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize(1); // no 'shadow'

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.90);

  gStyle->SetPalette(1); // rainbow colors

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();
 
  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //                topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );

  c1.Print( "thr563.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  gPad->Update();// required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read data:

  string HASH( "#" );

  cout << "try to open thr563.dat";
  ifstream D563( "thr563.dat" );
  if( !D563 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vthr;
  vector <double> vrow;
  vector <double> vres;

  while( D563.good() && ! D563.eof() ) {

    string rl;
    getline( D563, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int thr;
    double row;
    double res;

    strm >> thr;
    strm >> row;
    strm >> res;

    vthr.push_back(thr);
    vrow.push_back(row);
    vres.push_back(res);

    cout << vthr.size() << ". " << thr << "  " << row << endl;

  } // while lines

  int n563 = vthr.size();
  cout << n563 << " thr for chip 563" << endl;

  // subtract telescope:

  vector <double> vsig(n563);
  for( int i = 0; i < n563; ++i ) {
    vthr[i] = 0.7 + 0.55*(vthr[i]-1);
    vsig[i] = sqrt( vres[i]*vres[i] - 3.6*3.6 ); // sixdxc/2
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // row vs thr

  TH1F * hr = new
    TH1F( "hr",
	  "rows vs thr;threshold [ke];Lin cluster size [rows]",
	  15, 0, 5 ); // axis range
  hr->SetStats(kFALSE); // no statistics
  hr->SetMinimum(1);
  hr->SetMaximum(3);
  hr->SetNdivisions(-504,"y");
  hr->Draw();

  TGraph * gr = new TGraph( n563, &vthr[0], &vrow[0] );
  gr->SetLineColor(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1.5);
  gr->Draw("Pc"); // without axis option: overlay

  c1.Print( "thr563.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // sig vs thr

  TH1F * hs = new
    TH1F( "hs",
	  "resolution vs thr;threshold [ke];Lin resolution [#mum]",
	  15, 0, 5 ); // axis range
  hs->SetStats(kFALSE); // no statistics
  hs->SetMinimum(0);
  hs->SetMaximum(7);
  hs->Draw();

  TGraph * gs = new TGraph( n563, &vthr[0], &vsig[0] );
  gs->SetLineColor(2);
  gs->SetMarkerColor(2);
  gs->SetMarkerStyle(20);
  gs->SetMarkerSize(1.5);
  gs->Draw("Pc"); // without axis option: overlay

  c1.Print( "thr563.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "thr563.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf thr563.ps" );
  ierr = system( "rm -f  thr563.ps" );
  cout << "evince thr563.pdf" << endl;

}
