
// Daniel Pitzl, Dec 2018
// 509 Lin Q0 vs R
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
  gStyle->SetTitleOffset( 1.7, "y" );

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

  c1.Print( "q0vsr.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.83, 0.3, 0.97, 0.89 );

  // 800 V:

  TFile * f800 = TFile::Open( "scopeRD35116.root" );
  linqvsr->SetStats(0);
  linqvsr->SetTitle( "RD53 Lin 8#upoint10^{15 }p/cm^{2} " );
  linqvsr->GetXaxis()->SetTitle( "distance from irradiation center [mm]" );
  linqvsr->SetLineColor(1);
  linqvsr->SetMinimum(0);
  linqvsr->SetMaximum(12);
  linqvsr->Draw();
  lgnd->AddEntry( linqvsr, " 800 V", "l" );

  // 600 V:

  TFile * f600 = TFile::Open( "scopeRD35121.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(2);
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 600 V", "l" );

  // 500 V:

  TFile * f500 = TFile::Open( "scopeRD35122.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(4);
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 500 V", "l" );

  // 400 V:

  TFile * f400 = TFile::Open( "scopeRD35123.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(8);
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 400 V", "l" );

  // 300 V:

  TFile * f300 = TFile::Open( "scopeRD35125.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(6);
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 300 V", "l" );

  // 250 V:

  TFile * f250 = TFile::Open( "scopeRD35126.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(7); // green
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 250 V", "l" );

  // 200 V:

  TFile * f200 = TFile::Open( "scopeRD35127.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(94); // orange
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 200 V", "l" );

  // 150 V:

  TFile * f150 = TFile::Open( "scopeRD35131.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(28); // braun
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 150 V", "l" );

  // 100 V:

  TFile * f100 = TFile::Open( "scopeRD35132.root" );
  linqvsr->SetStats(0);
  linqvsr->SetLineColor(9); // blue
  linqvsr->Draw( "same" );
  lgnd->AddEntry( linqvsr, " 100 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "q0vsr.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "q0vsr.pdf]" ); // ] closes file
  cout << "evince q0vsr.pdf" << endl;

}
