%% Question 2
clear

%% Part a

ep =1; % initialising epsilon to 1

for n=1:100  % perfom loop
    n;
    ep(n+1) =  (ep(n)/2);
    if ((ep(n+1)/2) + 1) == 1
        break
    
    end

end
epa = ep(end); % print the value of epsilon that terminates the loop

stringa = sprintf('%.22f',epa); % converts into 22 digit string & prints in workspace.
fprintf('%.22f',epa); % prints the epsilon in the command window as a string

% n is printed in the workspace too.

%% Part b
 epb = 1;
 
while (epb/2) + 1 > 1
    epb = epb/2;  
end

stringb = sprintf('%.22f',epb);
fprintf('%.22f',epb);
%% Part c

epc = 1/2^(n); % keeps as number, same as 2^-n

stringc = sprintf('%.22f',epc); % converts to 22 digit string
fprintf('%.22f',epc);

%% Part d

ep1 = sum(epc +1); % equal to double, 1.000
ep2 = sum((epc/2) +1); % equal to 1, integer

isequal(epa,epb,epc,eps) % test to make sure its all correct, 1 = true.
