
// Daniel Pitzl, Jun 2019
// scattering Al
// root -l scat13.C

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>

//------------------------------------------------------------------------------
//void scat13()
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
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );

  c1.Print( "scat13.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  gPad->Update();// required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read data:

  string HASH( "#" );

  cout << "try to open scat13.dat";
  ifstream D13( "scat13.dat" );
  if( !D13 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vd;
  vector <double> v12;
  vector <double> v24;
  vector <double> ve;

  while( D13.good() && ! D13.eof() ) {

    string rl;
    getline( D13, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double d;
    double t24;
    double t12;

    strm >> d;
    strm >> t24;
    strm >> t12;

    vd.push_back(d);
    v24.push_back(t24);
    v12.push_back(t12);
    ve.push_back(0.001*t12); // error

  } // while lines

  int n13 = vd.size();
  cout << "read " << n13 << " data lines" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // plot:

  TH1F * ht = new
    TH1F( "ht",
	  "scattering in Al;Al thickness [mm];RMS scattering angle [mrad]",
	  8, 0, 8 ); // axis range
  ht->SetStats(kFALSE); // no statistics
  ht->SetMinimum(0);
  ht->SetMaximum(4);
  ht->SetNdivisions( -504, "y" );
  ht->Draw();

  TGraphErrors * g24 = new TGraphErrors( n13, &vd[0], &v24[0], 0, &ve[0] );
  g24->SetMarkerColor(4);
  g24->SetMarkerStyle(20);
  g24->SetMarkerSize(1.2);

  TF1 * H24 = new
    TF1( "H24", "sqrt([0]*[0]+pow(13.6/2.27*sqrt(x/[1])*(1+0.038*log((x+0.001)/[1])),2))", 0, 8 );
  H24->SetNpx(500);
  H24->SetLineWidth(2);
  H24->SetParameter( 0, 0.2 );
  H24->SetParameter( 1, 89 ); // X0 [mm]
  H24->SetLineColor(kCyan);
  //H24->Draw("same");

  cout << "Highland fit 2.4 GeV:" << endl;
  g24->Fit( "H24", "r" );
  g24->Draw("P"); // without axis option: overlay

  TGraphErrors * g12 = new TGraphErrors( n13, &vd[0], &v12[0], 0, &ve[0] );
  g12->SetMarkerColor(2);
  g12->SetMarkerStyle(21);
  g12->SetMarkerSize(1.2);

  TF1 * H12 = new
    TF1( "H12", "sqrt([0]*[0]+pow(13.6/1.2*sqrt(x/[1])*(1+0.038*log((x+0.001)/[1])),2))", 0, 8 );
  H12->SetNpx(500);
  H12->SetLineWidth(2);
  H12->SetParameter( 0, 0.4 );
  H12->SetParameter( 1, 89 ); // X0 [mm]
  H12->SetLineColor(kMagenta);
  //H12->Draw("same");

  cout << "Highland fit 1.2 GeV:" << endl;
  g12->Fit( "H12", "r" );
  g12->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.70, 0.17, 0.93, 0.32 );
  lgnd->AddEntry( g12, "1.2 GeV", "p" );
  lgnd->AddEntry( g24, "2.4 GeV", "p" );
  lgnd->Draw( "same" );

  c1.Print( "scat13.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "scat13.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf scat13.ps" );
  ierr = system( "rm -f  scat13.ps" );
  cout << "evince scat13.pdf" << endl;

}
