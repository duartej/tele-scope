
// Daniel Pitzl, Nov 2018
// RFD53A turn scan
// root -l turn524.C

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>

//------------------------------------------------------------------------------
void turn524()
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

  gStyle->SetTitleOffset( 1.5, "x" );
  gStyle->SetTitleOffset( 1.5, "y" );

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
  //                topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );

  c1.Print( "turn524.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read data:

  string HASH( "#" );

  cout << "try to open turn524.dat";
  ifstream D524( "turn524.dat" );
  if( !D524 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> turn;
  vector <double> linq;
  vector <double> difq;
  vector <double> vnpx;
  vector <double> vrxc;

  while( D524.good() && ! D524.eof() ) {

    string rl;
    getline( D524, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int run;
    double trn;
    double lin;
    double dif;
    double npx;
    double dxc;

    strm >> run;
    strm >> trn;
    strm >> lin;
    strm >> dif;
    strm >> npx;
    strm >> dxc;

    turn.push_back(trn);
    linq.push_back(lin);
    difq.push_back(dif);
    vnpx.push_back(npx);
    vrxc.push_back(dxc);

    cout << turn.size() << ". " << trn << "  " << lin << "  " << dif << endl;

  } // while lines

  int n524 = turn.size();
  cout << n524 << " turns for chip 524" << endl;

  vector <double> vsig(n524);
  for( int i = 0; i < n524; ++i ) {
    vsig[i] = sqrt( vrxc[i]*vrxc[i] - 6.5*6.5 ); // subtract track sixdxc/2
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // plot:

  TH1F * hq = new
    TH1F( "hq",
	  "cluster ToT vs turn;turn angle [deg];<cluster charge> [ToT]",
	  36, -6, 30 ); // axis range
  hq->SetStats(kFALSE); // no statistics
  hq->SetMinimum(10);
  hq->SetMaximum(25);
  hq->Draw();

  // create a TF1 for fit:

  TF1 * PathLength = new TF1( "Fpath", "[0]/cos(0.017453*x)", -5, 29 ); // deg->rad, range

  PathLength->SetNpx(500);
  PathLength->SetLineWidth(2);

  PathLength->SetParameter( 0, linq[18] ); // Q0
  PathLength->SetLineColor(kCyan);
  PathLength->Draw( "same" );
  c1.Update();
  /*
  PathLength->SetParameter( 0, difq[18] ); // Q0
  PathLength->SetLineColor(kMagenta);
  PathLength->Draw( "same" );
  c1.Update();
  */
  // Lin:

  TGraph * gql = new TGraph( n524, &turn[0], &linq[0] );
  gql->SetMarkerColor(4);
  gql->SetMarkerStyle(20);
  gql->SetMarkerSize(1.5);
  gql->Draw("P"); // without axis option: overlay

  // Diff:

  TGraph * gqd = new TGraph( n524, &turn[0], &difq[0] );
  gqd->SetMarkerColor(2);
  gqd->SetMarkerStyle(21);
  gqd->SetMarkerSize(1.5);
  gqd->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.70, 0.17, 0.93, 0.32 );
  lgnd->AddEntry( gqd, "524 Dif ", "p" );
  lgnd->AddEntry( gql, "524 Lin ", "p" );
  lgnd->Draw( "same" );

  c1.Print( "turn524.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // npx vs turn

  TH1F * hn = new
    TH1F( "hn",
	  "cluster size vs turn;turn angle [deg];Lin cluster size [pixels]",
	  36, -6, 30 ); // axis range
  hn->SetStats(kFALSE); // no statistics
  hn->SetMinimum(1);
  hn->SetMaximum(3);
  hn->Draw();

  TGraph * gn = new TGraph( n524, &turn[0], &vnpx[0] );
  gn->SetMarkerColor(1);
  gn->SetMarkerStyle(21);
  gn->SetMarkerSize(1.5);
  gn->Draw("P"); // without axis option: overlay

  c1.Print( "turn524.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // sig vs turn

  TH1F * hs = new
    TH1F( "hs",
	  "resolution vs turn;turn angle [deg];Lin resolution [#mum]",
	  36, -6, 30 ); // axis range
  hs->SetStats(kFALSE); // no statistics
  hs->SetMinimum(0);
  hs->SetMaximum(15);
  hs->Draw();

  TGraph * gs = new TGraph( n524, &turn[0], &vsig[0] );
  gs->SetMarkerColor(2);
  gs->SetMarkerStyle(22);
  gs->SetMarkerSize(1.5);
  gs->Draw("P"); // without axis option: overlay

  c1.Print( "turn524.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "turn524.pdf]" ); // ] closes file
  cout << "evince turn524.pdf" << endl;

}
