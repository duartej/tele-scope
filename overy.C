
// Feb 2018
// scopes.cc: overlay triy for alignment check
// rlq scopes34500.root
// .x overy.C
{
  gStyle->SetOptStat(111);
  //gStyle->SetStatY(0.90);
  //gStyle->SetStatW(0.30);
  //gStyle->SetStatH(0.30);

  triyc->Draw();
  N0 = triyc->Integral();
  gPad->SetGridx();

  triyclk->Draw("esame");
  N1 = triyclk->Integral();
  cout << double(N1)/N0 << endl;
}
