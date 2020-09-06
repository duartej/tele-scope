
// Daniel Pitzl, May 2019
// 511i eff vs xm bias scan
// root -l effvsxmbias.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.4, "y" );

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
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "effvsxmbias.pdf[", "pdf" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.78, 0.17, 0.93, 0.48 );

  // 800 V

  TFile * f33 = TFile::Open( "scopeRD35636.root" );
  effvsxm->SetStats(0);
  effvsxm->SetTitle( "511 Lin 7#upoint10^{15 }n_{eq}/cm^{2} " );
  effvsxm->SetLineColor(1);
  effvsxm->SetMinimum(0.6);
  effvsxm->SetMaximum(1.0);
  effvsxm->SetNdivisions( -1004, "y" );
  effvsxm->Draw();
  lgnd->AddEntry( effvsxm, " 754 V", "l" );

  // bias 700

  TFile * f27 = TFile::Open( "scopeRD35638.root" );
  effvsxm->SetStats(0);
  effvsxm->SetLineColor(2);
  effvsxm->Draw( "same" );
  lgnd->AddEntry( effvsxm, " 662 V", "l" );

  // bias 600

  TFile * f22 = TFile::Open( "scopeRD35639.root" );
  effvsxm->SetStats(0);
  effvsxm->SetLineColor(15); // grey
  effvsxm->SetLineColor(51); // violet
  effvsxm->Draw( "same" );
  lgnd->AddEntry( effvsxm, " 570 V", "l" );

  // bias 500

  TFile * f17 = TFile::Open( "scopeRD35640.root" );
  effvsxm->SetStats(0);
  effvsxm->SetLineColor(92); // gold
  effvsxm->SetLineColor(95); // orange
  effvsxm->Draw( "same" );
  lgnd->AddEntry( effvsxm, " 475 V", "l" );

  // bias 400

  TFile * f12 = TFile::Open( "scopeRD35641.root" );
  effvsxm->SetStats(0);
  effvsxm->SetLineColor(8);
  effvsxm->SetLineColor(416+3); // dark green
  effvsxm->Draw( "same" );
  lgnd->AddEntry( effvsxm, " 381 V", "l" );

  // bias 300

  TFile * f07 = TFile::Open( "scopeRD35642.root" );
  effvsxm->SetStats(0);
  effvsxm->SetLineColor(4);
  effvsxm->Draw( "same" );
  lgnd->AddEntry( effvsxm, " 286 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "effvsxmbias.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "effvsxmbias.pdf]" ); // ] closes file
  cout << "evince effvsxmbias.pdf" << endl;

}
