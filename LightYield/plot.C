
{

    TGraph * testcan = new TGraph("testcan.dat"); 
    TGraph * prototype = new TGraph("prototype.dat");
    TF1 * func = new TF1("func","[0]*x - [1]*(1-TMath::Exp(-[2]*TMath::Power(x,[3])))",0,20);

    //func->SetParameters(0.55322,1.58611,0.29188,0.98036);
    func->SetParameters(0.67353 , 2.64302 , 0.21628 , 1.00408); // with offset0 test5 parameters

    testcan->SetMarkerStyle(25);
    prototype->SetMarkerStyle(26);

    testcan->SetMarkerSize(0.8);
    prototype->SetMarkerSize(0.8);
    
    testcan->SetMarkerColor(kBlack);
    prototype->SetMarkerColor(kBlue);
    func->SetLineColor(kRed);

    func->Draw();
    prototype->Draw("p same");
    testcan->Draw("p same");

}
