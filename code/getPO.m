% [coeffs, lbls] = getPO(rho, threshold);
%
% getPO takes a density matrix and decomposes it in terms
% of the product operators.  The function returns only
% those PO terms with coefficients greater than threshold.
%
% returns the coefficients and the coefficient labels...
function  [coeffs, lbls] = getPO(rho, threshold)

I{1}=.5*[1,0;0,1]; 
I{2}=.5*[0,1;1,0]; 
I{3}=.5*[0,-1i;1i,0]; 
I{4}=.5*[1,0;0,-1];

text{1}='I';
text{2}='X';
text{3}='Y';
text{4}='Z';

index = ones(1,log2(length(rho)));

%nn indexes coeffecients that are returned.
nn = 1; 

% Create PO Basis and calculates projection of rho
% onto each of the basis states.
for j=1:length(rho)^2,
  
  if index(1) == -10,
    disp('ERROR: index reached maximum value and still looping')
    break;
  end
  
  PO = 1;
  POText = '';
  
  for k=1:length(index),
    PO = kron(PO,I{index(k)});
    POText = strcat(POText,text{index(k)});
  end
  
  POCoeffs(j) = trace(rho*PO);
  
  if abs(POCoeffs(j)) > threshold,
    coeffs(nn) = POCoeffs(j);
    lbls(nn,:) = POText;
    disp(sprintf('%s:\t%1.4f+%1.4fi\n',lbls(nn,:),real(coeffs(nn)),imag(coeffs(nn))));
    nn=nn+1;

  end

  index = incIndex(index);

end

return;


%%%%
function index=incIndex(index)

lenIndex = length(index);

index(lenIndex)=index(lenIndex)+1;

if index(lenIndex) > 4 && lenIndex > 1, 
  index(lenIndex) = 1;
  index(1:lenIndex-1) = incIndex(index(1:lenIndex-1));
elseif lenIndex == 1 && index(lenIndex) > 4 , 
  index = -10;
end











