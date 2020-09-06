
// Daniel Pitzl, Mar 2020
// RD53A vs fluence
// root -l v99.C

//------------------------------------------------------------------------------
//void flu()
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleOffset( 1.4, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLineWidth(1);// frames
  //gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineColor(1); // 1=black
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
  //             topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );

  c1.Print( "v99.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // proton irradiated no bias dot:

  const Int_t np = 9; // order by facility for color code
  //                  PS   Ka   Ka   Ka   Ka   Ka   Bh   Bh   Bh
  Float_t scm[np] = { 509, 547, 578, 542, 577, 589, 564, 565, 563 };
  Float_t flu[np] = { 5.2, 4.4, 5.0, 7.4, 9.9,21.4, 5.0, 5.0,10.1 };
  Float_t w99[np] = { 380, 330, 220, 400, 440, 660, 420, 405, 460 };
  Float_t v99[np] = {-242, 330,-220, 380, 440, 660, 390, 430, 450 };

  Float_t efl[np]; // dosimetry error
  Float_t ev9[np]; // v99 error
  for( int i = 0; i < np; ++i ) {
    efl[i] = 0.07*flu[i];
    ev9[i] = 10; // [V]
  }

  // multiple bias scans:

  const int n09{2}; // PS
  Float_t f09[n09] { 5.2, 5.2 };
  Float_t v09[n09] { 280, 380 };

  const int n47{1}; // Ka
  Float_t f47[n47] { 4.4 };
  Float_t v47[n47] { 330 };

  const int n42{2}; // Ka
  Float_t f42[n42] { 7.4, 7.4 };
  Float_t v42[n42] { 380, 400 };

  const int n77{1}; // Ka spark
  Float_t f77[n77] { 10.1 };
  Float_t v77[n77] { 440 };

  const int n89{1}; // Ka
  Float_t f89[n89] { 21.4 };
  Float_t v89[n89] { 660 };

  const int n64{3}; // Bh
  Float_t f64[n64] { 4.9, 4.9, 4.9 };
  Float_t v64[n64] { 390, 450, 535 };

  const int n65{3}; // Bh
  Float_t f65[n65] { 5.1, 5.1, 5.1 };
  Float_t v65[n65] { 430, 490, 540 };

  const int n63{4}; // Bh
  Float_t f63[n63] { 9.9, 9.9, 9.9, 9.9 };
  Float_t v63[n63] { 450, 560, 665, 790 };

  // neutrons:

  // 194 = FDB150P R4S 100x25 P1 p-stop default rot90 with 4E15 n/cm2
  // 202 = FDB150Y R4S 50x50 Y2 p-spray default on straight  4E15 n  no bias scan done

  // 196 = FDB150P R4S 100x25 P1 p-stop default rot90 after 8E15 n/cm2
  // 206 = FDB150Y R4S 100x25 Y3 p-spray max imp rot90  8E15 n
  // 198 = FDB150P R4S 50x50 P1 p-stop default straight  8E15 n

  // 197 = FDB150P R4S 100x25 P7 p-stop max imp rot90 with 1.6E16 n/cm2
  // 207 = FDB150Y R4S 100x25 Y2 p-spray default rot90  1.6E16 n
  // 201 = FDB150P R4S 100x25 P1 p-stop default on straight PCB with 1.6E16 n

  const Int_t nn = 7;
  Float_t scn[nn] = {  194, 196, 206, 198,  197,  207,  201 }; // R4S single chip modules
  Float_t fln[nn] = {  3.6, 7.2, 7.2, 7.2, 14.4, 14.4, 14.4 };
  Float_t v9n[nn] = {  100, 240, 160, 270,  480,  410,  550 };

  Float_t efn[nn]; // dosimetry error
  Float_t e9n[nn]; // v99 error
  for( int i = 0; i < nn; ++i ) {
    efn[i] = 0.07*fln[i];
    e9n[i] = 10; // [V]
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // multiple v99 vs fluence

  hm = new
    TH1F( "hm",
	  "increasing thresholds;fluence [10^{15}n_{eq}/cm^{2}];bias for 99% efficiency [V]",
	  23, 0, 23 ); // axis range

  hm->SetStats(kFALSE); // no statistics
  hm->SetMinimum(0);
  hm->SetMaximum(800);
  hm->Draw();

  TLegend * lm = new TLegend( 0.76, 0.16, 0.94, 0.60 );

  g09 = new TGraph( n09, f09, v09 );
  g09->SetLineWidth(2);
  g09->SetLineColor(4);
  g09->SetMarkerColor(4);
  g09->SetMarkerStyle(20);
  g09->SetMarkerSize(1.5);
  g09->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g09, "509 PS", "pl" );

  g47 = new TGraph( n47, f47, v47 );
  g47->SetLineWidth(2);
  g47->SetLineColor(kGreen+1);
  g47->SetMarkerColor(kGreen+1);
  g47->SetMarkerStyle(20);
  g47->SetMarkerSize(1.5);
  g47->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g47, "547 Ka", "pl" );

  g42 = new TGraph( n42, f42, v42 );
  g42->SetLineWidth(2);
  g42->SetLineColor(kGreen+2);
  g42->SetMarkerColor(kGreen+2);
  g42->SetMarkerStyle(20);
  g42->SetMarkerSize(1.5);
  g42->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g42, "542 Ka", "pl" );

  g77 = new TGraph( n77, f77, v77 );
  g77->SetLineWidth(2);
  g77->SetLineColor(kGreen+3);
  g77->SetMarkerColor(kGreen+3);
  g77->SetMarkerStyle(20);
  g77->SetMarkerSize(1.5);
  g77->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g77, "577 Ka", "pl" );

  g89 = new TGraph( n89, f89, v89 );
  g89->SetLineWidth(2);
  g89->SetLineColor(kGreen+4);
  g89->SetMarkerColor(kGreen+4);
  g89->SetMarkerStyle(20);
  g89->SetMarkerSize(1.5);
  g89->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g89, "589 Ka", "pl" );

  g64 = new TGraph( n64, f64, v64 );
  g64->SetLineWidth(2);
  g64->SetLineColor(618); // magenta
  g64->SetMarkerColor(618);
  g64->SetMarkerStyle(20);
  g64->SetMarkerSize(1.5);
  g64->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g64, "564 Bh", "pl" );

  g65 = new TGraph( n65, f65, v65 );
  g65->SetLineWidth(2);
  g65->SetLineColor(2);
  g65->SetMarkerColor(2);
  g65->SetMarkerStyle(20);
  g65->SetMarkerSize(1.5);
  g65->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g65, "565 Bh", "pl" );

  g63 = new TGraph( n63, f63, v63 );
  g63->SetLineWidth(2);
  g63->SetLineColor(634); // dark red
  g63->SetMarkerColor(634);
  g63->SetMarkerStyle(20);
  g63->SetMarkerSize(1.5);
  g63->Draw("Pc"); // without axis option: overlay
  lm->AddEntry( g63, "563 Bh", "pl" );

  lm->Draw( "same" );

  c1.Print( "v99.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // v99 vs fluence

  h99 = new
    TH1F( "h99",
	  //"CMS pixel performance;fluence [10^{15}n_{eq}/cm^{2}];bias for 99% efficiency [V]",
	  ";fluence [10^{15}n_{eq}/cm^{2}];bias for 99% efficiency [V]",
	  24, 0, 24 ); // axis range

  h99->SetStats(kFALSE); // no statistics
  h99->SetMinimum(0);
  h99->SetMaximum(800);
  h99->Draw();

  TLegend * lgnde = new TLegend( 0.58, 0.16, 0.94, 0.36 );

  gPS = new TGraphErrors( 1, &flu[0], &v99[0], &efl[0], &ev9[0] );
  gPS->SetLineWidth(2);
  gPS->SetLineColor(4);
  gPS->SetMarkerColor(4);
  gPS->SetMarkerStyle(20);
  gPS->SetMarkerSize(1.5);
  gPS->Draw("P"); // without axis option: overlay
  lgnde->AddEntry( gPS, "RD53A PS ", "pl" );

  gKa = new TGraphErrors( 5, &flu[1], &v99[1], &efl[1], &ev9[1] );
  gKa->SetLineWidth(2);
  gKa->SetLineColor(kGreen+3);
  gKa->SetMarkerColor(kGreen+3);
  gKa->SetMarkerStyle(23);
  gKa->SetMarkerSize(1.5);
  gKa->Draw("P"); // without axis option: overlay
  lgnde->AddEntry( gKa, "RD53A Ka ", "pl" );

  gBh = new TGraphErrors( 3, &flu[6], &v99[6], &efl[6], &ev9[6] );
  gBh->SetLineWidth(2);
  gBh->SetLineColor(2);
  gBh->SetMarkerColor(2);
  gBh->SetMarkerStyle(22);
  gBh->SetMarkerSize(1.5);
  gBh->Draw("P"); // without axis option: overlay
  lgnde->AddEntry( gBh, "RD53A Bh ", "pl" );
  /*
  gn = new TGraphErrors( nn, fln, v9n, efn, e9n );
  gn->SetLineWidth(2);
  gn->SetLineColor(kOrange-3);
  gn->SetMarkerColor(kOrange-3);
  gn->SetMarkerStyle(20);
  gn->SetMarkerSize(1.5);
  gn->Draw("P"); // without axis option: overlay
  lgnde->AddEntry( gn , "R4S n  ", "pl" );
  */
  lgnde->Draw( "same" );

  c1.Print( "v99.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "v99.pdf]" ); // ] closes file
  cout << "evince v99.pdf" << endl;

}
