
// Daniel Pitzl, Jun 2019
// scattering CFK
// root -l scat6.C

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>

//------------------------------------------------------------------------------
//void scat6()
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );

  c1.Print( "scat6.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  gPad->Update();// required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read data:

  string HASH( "#" );

  cout << "try to open scat6.dat";
  ifstream D6( "scat6.dat" );
  if( !D6 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vd;
  vector <double> v12;
  vector <double> v08;
  vector <double> ve;

  while( D6.good() && ! D6.eof() ) {

    string rl;
    getline( D6, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double d;
    double t08;
    double t12;

    strm >> d;
    strm >> t12;
    strm >> t08;

    vd.push_back(d);
    v08.push_back(t08);
    v12.push_back(t12);
    ve.push_back(0.001*t12); // error

  } // while lines

  int n6 = vd.size();
  cout << "read " << n6 << " data lines" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // plot:

  TH1F * ht = new
    TH1F( "ht",
	  "scattering in CFK;CFK thickness [mm];RMS scattering angle [mrad]",
	  3, 0, 3 ); // axis range
  ht->SetStats(kFALSE); // no statistics
  ht->SetMinimum(0);
  ht->SetMaximum(1.6);
  //ht->SetNdivisions( -504, "y" );
  ht->Draw();

  TGraphErrors * g08 = new TGraphErrors( n6, &vd[0], &v08[0], 0, &ve[0] );
  g08->SetMarkerColor(4);
  g08->SetMarkerStyle(20);
  g08->SetMarkerSize(1.2);

  TF1 * H08 = new
    TF1( "H08", "sqrt([0]*[0]+pow(13.6/0.84*sqrt(x/[1])*(1+0.038*log((x+0.001)/[1])),2))", 0, 3 );
  H08->SetNpx(500);
  H08->SetLineWidth(2);
  H08->SetParameter( 0, 0.5 );
  H08->SetParameter( 1, 360 ); // X0 [mm]
  H08->SetLineColor(kCyan);
  //H08->Draw("same");

  cout << "Highland fit 0.8 GeV:" << endl;
  g08->Fit( "H08", "r" );
  g08->Draw("P"); // without axis option: overlay

  TGraphErrors * g12 = new TGraphErrors( n6, &vd[0], &v12[0], 0, &ve[0] );
  g12->SetMarkerColor(2);
  g12->SetMarkerStyle(21);
  g12->SetMarkerSize(1.2);

  TF1 * H12 = new
    TF1( "H12", "sqrt([0]*[0]+pow(13.6/1.2*sqrt(x/[1])*(1+0.038*log((x+0.001)/[1])),2))", 0, 3 );
  H12->SetNpx(500);
  H12->SetLineWidth(2);
  H12->SetParameter( 0, 0.4 );
  H12->SetParameter( 1, 360 ); // X0 [mm]
  H12->SetLineColor(kMagenta);
  //H12->Draw("same");

  cout << "Highland fit 1.2 GeV:" << endl;
  g12->Fit( "H12", "r" );
  g12->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.70, 0.17, 0.93, 0.32 );
  lgnd->AddEntry( g08, "0.8 GeV", "p" );
  lgnd->AddEntry( g12, "1.2 GeV", "p" );
  lgnd->Draw( "same" );

  c1.Print( "scat6.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "scat6.pdf]" ); // ] closes file
  cout << "evince scat6.pdf" << endl;

}
