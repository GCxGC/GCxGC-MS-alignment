function [MSvalueboxRounded,MSintboxRounded] = MSdataRound_v2(MSvaluebox,...
    MSintbox,Precision)

tic

% Roounds MS m/z data according to precision given in Precision (e.g.
% 0.001). 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% See license terms stated in LICENSE.txt
% Authors : Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

 % Round m/z values:
  MSvalueboxII= round(MSvaluebox*(1/Precision))*Precision;
MSintboxII = MSintbox;
 
 % We will now just , for each pixel, only keep once each m/z value, and sum
% up the corresponding intensity values in case the m/z value was present
% more than once:
% A & B will contain the data, without multiple ion occurences:
A = zeros(size(MSvalueboxII));
B = zeros(size(MSvalueboxII));
% Let's go through each line (kt is the line number):
for kt=1:size(MSvalueboxII,1)
%     First element is just the first element in the aligned MS data:
    B(kt,1) = MSvalueboxII(kt,1);
    A(kt,1) = MSintboxII(kt,1);
%     "Cnt" will refer to the elements in the current line of the new 
% A & B matrices:
    Cnt = 1;
%     Let's now go through the current line (so go through it column by 
% column):
    for rr=1:size(MSvalueboxII,2)
%         If we already have the MS ion into B, we just add the intensity
%         value at the corresponding place in A:
        if MSvalueboxII(kt,rr) == B(kt,Cnt)
            A(kt,Cnt) = A(kt,Cnt) + MSintboxII(kt,rr);
%             Else, we just keep the values at the next position in A & B.
        else
            Cnt = Cnt+1;
            B(kt,Cnt) = MSvalueboxII(kt,rr);
            A(kt,Cnt) = MSintboxII(kt,rr);
        end
        
    end
end

% We will again remove useless zeros:
MaxNotZero2=max(sum(B~=0,2));
% MaxNotZero2=max(sum(MSvalueboxII~=0,2))

% Remove the useless zeros from the result:
B = B(:,1:MaxNotZero2);
A = A(:,1:MaxNotZero2);
MSvalueboxRounded=B;
MSintboxRounded=A;
% MSvalueboxRounded=MSvalueboxII(:,1:MaxNotZero2);
% MSintboxRounded=MSintboxII(:,1:MaxNotZero2);

disp('Rounding of m/z data achieved in:')
toc


