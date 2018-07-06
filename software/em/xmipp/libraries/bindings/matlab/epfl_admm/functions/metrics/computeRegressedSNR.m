function SNR = computeRegressedSNR(rec,oracle)

sumP    =        sum(sum(oracle))                  ;
sumI    =        sum(sum(rec))                      ;
sumIP   =        sum(sum(oracle.*rec))              ;
sumI2   =        sum(sum(rec.^2))                   ;
A       =        [sumI2,sumI;sumI,size(oracle,1)*size(oracle,2)]        ;
b       =        [sumIP;sumP]                      ;
c       =        (A)\b                             ;
rec     =        c(1)*rec+c(2)                     ;
err     =        sum(sum((oracle-rec).^2))         ;
SNR     =        10*log10(sum(sum(oracle.^2))/err) ;