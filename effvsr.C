
// Daniel Pitzl, Dec 2018
// 509 Lin eff vs R
// root -l effvsr.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //               topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "effvsr.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.83, 0.3, 0.97, 0.89 );

  // 800 V:

  TFile * f800 = TFile::Open( "scopeRD35116.root" );
  effvsr->SetStats(0);
  effvsr->SetTitle( "RD53 Lin 8#upoint10^{15 }p/cm^{2} " );
  effvsr->GetXaxis()->SetTitle( "distance from irradiation center [mm]" );
  effvsr->SetLineColor(1);
  effvsr->SetMinimum(0.5);
  effvsr->SetMaximum(1.0);
  effvsr->Draw();
  lgnd->AddEntry( effvsr, " 800 V", "l" );

  // 600 V:

  TFile * f600 = TFile::Open( "scopeRD35121.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(2);
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 600 V", "l" );

  // 500 V:

  TFile * f500 = TFile::Open( "scopeRD35122.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(4);
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 500 V", "l" );

  // 400 V:

  TFile * f400 = TFile::Open( "scopeRD35123.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(8);
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 400 V", "l" );

  // 300 V:

  TFile * f300 = TFile::Open( "scopeRD35125.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(6);
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 300 V", "l" );

  // 250 V:

  TFile * f250 = TFile::Open( "scopeRD35126.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(7); // green
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 250 V", "l" );

  // 200 V:

  TFile * f200 = TFile::Open( "scopeRD35127.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(94); // orange
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 200 V", "l" );

  // 150 V:

  TFile * f150 = TFile::Open( "scopeRD35131.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(28); // braun
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 150 V", "l" );

  // 100 V:

  TFile * f100 = TFile::Open( "scopeRD35132.root" );
  effvsr->SetStats(0);
  effvsr->SetLineColor(9); // blue
  effvsr->Draw( "same" );
  lgnd->AddEntry( effvsr, " 100 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "effvsr.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "effvsr.pdf]" ); // ] closes file
  cout << "evince effvsr.pdf" << endl;

}
