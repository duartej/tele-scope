
// Daniel Pitzl, Mar 2020
// Datura triplet residuals vs p
// root -l tri.C

//------------------------------------------------------------------------------
//void tri()
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

  c1.Print( "tri.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read data:

  string HASH( "#" );

  cout << "try to open tri.dat";
  ifstream Dt( "tri.dat" );
  if( !Dt ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vp;
  vector <double> v20;
  vector <double> v74;
  vector <double> vep;
  vector <double> ver;

  while( Dt.good() && ! Dt.eof() ) {

    string rl;
    getline( Dt, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double p;
    double r20;
    double r74;

    strm >> p;
    strm >> r20;
    strm >> r74;

    vp.push_back(p);
    v20.push_back(r20);
    v74.push_back(r74);
    vep.push_back(0.15); // error
    ver.push_back(0.10); // error

  } // while lines

  int np = vp.size();
  cout << "read " << np << " data lines" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // plot:

  TH1F * hr = new
    TH1F( "hr",
	  "Datura triplet residual;beam momentum [GeV];triplet residual [#mum]",
	  6, 0, 6 ); // axis range
  hr->SetStats(kFALSE); // no statistics
  hr->SetMinimum( 0);
  hr->SetMaximum(11);
  //hr->SetNdivisions( -504, "y" );
  hr->Draw();

  TGraphErrors * g20 = new TGraphErrors( np, &vp[0], &v20[0], &vep[0], &ver[0] );
  g20->SetMarkerColor(4);
  g20->SetMarkerStyle(20);
  g20->SetMarkerSize(1.2);

  TF1 * f20 = new
    TF1( "f20", "sqrt([0]/x*[0]/x+1.5*[1]*[1])", 0.5, 6 );
  f20->SetNpx(500);
  f20->SetLineWidth(2);
  f20->SetParameter( 0, 5 );
  f20->SetParameter( 1, 3.6 );
  f20->SetLineColor(kCyan);
  //f20->Draw("same");

  cout << endl << "fit 20 mm:" << endl;
  g20->Fit( "f20", "r" );
  g20->Draw("P"); // without axis option: overlay

  TGraphErrors * g74 = new TGraphErrors( np, &vp[0], &v74[0], &vep[0], &ver[0] );
  g74->SetMarkerColor(2);
  g74->SetMarkerStyle(21);
  g74->SetMarkerSize(1.2);

  TF1 * f74 = new
    TF1( "f74", "sqrt([0]/x*[0]/x+1.5*[1]*[1])", 0.5, 6 );
  f74->SetNpx(500);
  f74->SetLineWidth(2);
  f74->SetParameter( 0, 11 );
  f74->SetParameter( 1, 3.5 );
  f74->SetLineColor(kMagenta);
  //f74->Draw("same");

  cout << endl << "fit 74 mm:" << endl;
  g74->Fit( "f74", "r" );
  g74->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.70, 0.17, 0.93, 0.32 );
  lgnd->AddEntry( g74, "74 mm", "p" );
  lgnd->AddEntry( g20, "20 mm", "p" );
  lgnd->Draw( "same" );

  c1.Print( "tri.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "tri.pdf]" ); // ] closes file
  cout << "evince tri.pdf" << endl;

}
