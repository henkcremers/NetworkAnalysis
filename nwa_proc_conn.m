function [Cproc,diagnostic] = nwa_proc_conn(C,varargin)
% ==================================================
% [Cproc] = nwa_reshape(C,varargin)
% EXAMPLES
% =================================================
nn = size(C);
Cproc = C;
diagnostic = [];

% nd = 1:nn(3);
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            
            % treshold with absolute values (a)
            case 'thresha'
                tval = varargin{i+1};
                loc = Cproc<tval;
                Cproc(loc) = 0;
                
            % treshold with relative value (r)
            case 'threshr'
                dat = nwa_reshape(Cproc,'mat2vec');
                dat = abs(dat);
                sortdat = sort(dat);
                thrP = varargin{i+1};
                thr = sortdat(ceil(length(sortdat)*(1-thrP)));
                %Cproc = zeros(size(C));
                loc = abs(C)>=thr;
                Cproc(~loc) = 0; 
                %C(loc);
                
            % binarize
            case 'bin'
                %binval = varargin{i+1};
                loc = Cproc~=0;
                Cproc(loc) = 1;
                %Cproc(~loc) = 1;
                
            % take absolute values
            case 'abs'
                Cproc = abs(Cproc);
                
            % make symmetric by averaging (e.g. for results of gppi)
            case 'avesym'
                Cflip = nwa_reshape(Cproc,'flipdiag');
                Ctemp(:,:,1) = Cproc;
                Ctemp(:,:,2) = Cflip;
                Cproc = mean(Ctemp,3);
                % cdiff = abs(Cproc-Cflip);% i think using Cproc is wrong,
                % this is already the average. it should be:
                cdiff = abs(Ctemp(:,:,1) - Ctemp(:,:,2)); 
                %cabs  = abs(Ctemp(:,:,1)) + abs(Ctemp(:,:,2));
                csum = Ctemp(:,:,1) + Ctemp(:,:,2);
                %c2mean = 2.*abs(Cproc);
                cdiag = abs(cdiff./csum);
                cdiag = nwa_reshape(cdiag,'mat2vec');

                % diagnostic = nanmean(cdiag);
                % only incorporate the "relevant" (highest 5%) values. 
                % THIS MAKES A BIG DIFFERENCE
                Cs = nwa_procmat(Cproc,'threshr',0.05,'bin');
                csvec = logical(nwa_reshape(Cs,'mat2vec'));
                diagnostic(1) = nanmean(cdiag(csvec));
                
                c1 = nwa_reshape(Ctemp(:,:,1),'mat2vec');
                c2 = nwa_reshape(Ctemp(:,:,2),'mat2vec');
                d = c1(csvec) - c2(csvec);
                diagnostic(2:3) =range(d);
                
           
            % set diagonal to zero    
            case 'diag0' 
                Cproc(logical(eye(nn(1)))) = 0;
             
             % scale values to a range of 0 - 1    
            case 'scale'
                Cproc = abs(Cproc);
                maxval = max(max(Cproc));
                loc = Cproc~=0;
                Cproc(loc) = Cproc(loc)./maxval;
                    
            otherwise
                error(['unknown input: ' arg])
                
                
        end
    end
end
return
