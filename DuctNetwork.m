classdef DuctNetwork < handle
    % DuctNetwork implements all algorithms for balancing duct network
    properties
        n; % number of nodes without ground
        n_NodeDescription; % cell array of size n by 1, text description on each node
        P; % array of size n by 1, pressure vector on each node , Pa
        b; % number of branches in the network
        b_FittingDescription; % cell array of size b by 1, text description on each branch
        Q; % array of size b by 1, flow rate vector on each branch, m^3/s
        A; % Associate matrix of size n by b
        % A(i,j)==1 means branch j leaves node i
        % A(i,j)==0 means branch j is not connected with node i
        % A(i,j)==-1 means branch j enters node i
        t; % dimension of null space
        U; % null space basis of matrix A, of size b by t, AU=0, U'U=1
        X; % internal state of duct system, of size t by 1
        X_lb; % lower boundary of possible internal state in simulation
        X_ub; % upper boundary of possible internal state in simulation
        b_Pdrop; % cell array of size b by 1 for the pressure drop in terms of q and s
        % if there exist multiple Pdrop functions, use a cell array to store
        % each of them only record the original Pdrop function for the
        % fitting, handle the fitting direction in pressure drop calculation
        b_dPdQ; % cell array of size b by 1 for the partial pressure drop over q in terms of q and s
        % only record the original dPdQ function for the fitting, handle the
        % fitting direction in pressure drop calculation
        b_dPdS; % cell array of size b by 1 for the partial pressure drop over s in terms of q and s
        % only record the original dPdS function for the fitting, handle the
        % fitting direction in pressure drop calculation
        b_Qidx; % cell array of size b by 1 for the index of dependent Q of Pdrop, dPdQ and dPdS functions, so use Q(abs(b_Qidx{b})).*sign(b_Qidx{b}) in these functions
        % Qidx{ii}[jj] can be negative, the sign represents the fitting
        % direction w.r.t. branch direction
        b_Sidx; % cell array of size b by 1 for the index of dependent S of Pdrop, dPdQ and dPdS functions, so use S(b_Sidx{b}) in these functions
        s; % number of parameters in the model
        s_ParamDescription; % cell array of size b by 1, text description on each parameter
        S; % array of size s by 1, parameter vector for whole system for identification
        S_lb; % lower boundary of parameter in simulation
        S_ub; % upper boundary of parameter in simulation
        S_scale; % scale of parameter for normalization in identification
        s_m; % array of size s by 1, s_m(i)=0 means no need to identify S(i), s_m(i)>0 is the multiplicity of this parameter that need to be identified
        s_MultiS; % cell array of size s by 1, containing values of multiplicities of each parameter in row vectors, Use [s_MultiS{:}]' to obtain all in one column vector
        cnfg; % number of configurations in identification
        cnfg_MultiSIdx; % cell array of size cnfg by 1, each is s by 1 colume vector of index of identified parameter in [s_MultiS{:}]'
        ob; % number of sensors used in each procedure
        ob_Uncertainty; % number array of sensor uncertainty, size ob by 1
        proc; % number of procedures in identification experiments
        proc_ConfigIdx; % vector to indicate the configuration ID used in each procedure, size proc by 1
        proc_ObPosMatrix; % cell array of Observation Position Matrix, size proc by 1, each is matrix of size ob by (n+b), branch first and then node
        proc_Measurements; % number array of Measurements record, size proc by ob
        d; % number of dampers
        d_Sidx; % array of parameters index for damper positions
        gdr; % number of terminals
        gdr_Bidx; % array of branch index for GDRs
        gdr_Nidx; % array of branch index for GDRs
        n_trail;
        UseNumericGrad; % the jacobian option for optimizer
        Display; % The display option for optimizer
        Max_Trails; % The maximum trail of iterations
        StateResidualControlScale; % The simulation control scale of states, multiplied on res_Sim
        
    end
    
    methods
        function obj = DuctNetwork()
            obj.n = 0;
            obj.n_NodeDescription = cell(0,1);
            obj.P = zeros(0,1);
            obj.b = 0;
            obj.b_FittingDescription = cell(0,1);
            obj.Q = zeros(0,1);
            obj.A = zeros(0,0);
            obj.t = 0;
            obj.U = zeros(0,0);
            obj.X = [];
            obj.X_ub = [];
            obj.X_lb = [];
            obj.b_Pdrop = cell(0,1);
            obj.b_dPdQ = cell(0,1);
            obj.b_dPdS = cell(0,1);
            obj.b_Qidx = cell(0,1);
            obj.b_Sidx = cell(0,1);
            obj.s = 0;
            obj.s_ParamDescription = cell(0,1);
            obj.S = zeros(0,1);
            obj.S_lb = zeros(0,1);
            obj.S_ub = zeros(0,1);
            obj.S_scale = ones(0,1);
            obj.s_m = zeros(0,1);
            obj.s_MultiS = cell(0,1);
            obj.cnfg=0;
            obj.cnfg_MultiSIdx = cell(0,1);
            obj.ob=0;
            obj.ob_Uncertainty = zeros(0,1);
            obj.proc=0;
            obj.proc_ConfigIdx = zeros(0,1);
            obj.proc_Measurements = zeros(0,0);
            obj.d = 0;
            obj.d_Sidx = [];
            obj.gdr = 0;
            obj.gdr_Bidx = [];
            obj.n_trail = 0;
            obj.UseNumericGrad = false;
            obj.Display = 'none';
            obj.Max_Trails = 10;
            obj.StateResidualControlScale = 1;
        end
        
        function Branch_Idx = AddBranch(obj, FromNode, ToNode, varargin)
            [~, b_Curr]=size(obj.A);
            Branch_Idx = b_Curr+1;
            SetNode(FromNode, 1);
            SetNode(ToNode, -1);
            [obj.n,obj.b]=size(obj.A);
            obj.U = 1e-1*null(obj.A);
            obj.t = size(obj.U,2);
            obj.X = ones(obj.t,1);
            obj.X_lb = -1e12*ones(obj.t,1);
            obj.X_ub = 1e12*ones(obj.t,1);
            obj.b_Pdrop{Branch_Idx,1}=cell(0,1);
            obj.b_dPdQ{Branch_Idx,1}=cell(0,1);
            obj.b_dPdS{Branch_Idx,1}=cell(0,1);
            obj.b_Qidx{Branch_Idx,1}=cell(0,1);
            obj.b_Sidx{Branch_Idx,1}=cell(0,1);
            obj.b_FittingDescription{Branch_Idx,1}=cell(0,1);
            if nargin>3
                obj.AddFitting(Branch_Idx,varargin{:});
            end
            function SetNode(Node, a)
                if ischar(Node), Node = find(cellfun(@(S)strcmp(S,Node), obj.n_NodeDescription),1); end
                if isempty(Node), Node=0; end
                if Node>0, obj.A(Node,Branch_Idx) = a; end
            end
        end
        
        function Node = AddNode(obj, varargin)
            Node = obj.n+1;
            obj.n = Node;
            if nargin>1
                obj.n_NodeDescription{Node,1} = varargin{1};
            end
        end
        
        function AddFitting(obj, Branch, Model, Param) % Branch has direction, positve means fitting direction is same as the branch direction, negative means opposite
            if isa(Model,'function_handle')
                if Model('Is_Junction')
                    % separate the junction into three/four junction branches, and add each of them separately, parameters are shared by all of them
                    CenterNode=find(arrayfun(@(i)isempty(setxor(find(obj.A(i,:)),abs(Branch))),1:obj.n),1);
                    if isempty(CenterNode)
                        disp('no such junction exist, fail to create junction'); return
                    end
                    Branch = -obj.A(CenterNode,abs(Branch)).*abs(Branch);
                    branch_function_handle = Model('Get_Branches');
                    branch_assignment = Model('Branch_Assignment');
                    param_assignment = Model('Parameter_Assignment');
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(ParamDescription));
                    for ii = 1:length(ParamDescription)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(strcmp(obj.s_ParamDescription,ParamDescription{ii}),1);
                            if isempty(s_Idx) % shared parameter is first included in the model
                                obj.s = obj.s+1; % add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; % add the parameter at the end
                            obj.s_ParamDescription{obj.s} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    for ii = 1:length(branch_function_handle)
                        BranchList = Branch(branch_assignment{ii});
                        ParamList = Param_Idx(param_assignment{ii});
                        PrimaryBranch = abs(BranchList(1));
                        obj.b_FittingDescription{PrimaryBranch}{end+1} = branch_function_handle{ii}('Model_Description');
                        obj.b_Pdrop{PrimaryBranch}{end+1} = @(q,s)branch_function_handle{ii}({'Pdrop','dPdQ','dPdS'},q,s);
                        obj.b_dPdQ{PrimaryBranch}{end+1} = @(q,s)branch_function_handle{ii}('dPdQ',q,s);
                        obj.b_dPdS{PrimaryBranch}{end+1} = @(q,s)branch_function_handle{ii}('dPdS',q,s);
                        obj.b_Qidx{PrimaryBranch}{end+1} = reshape(BranchList,[],1);
                        obj.b_Sidx{PrimaryBranch}{end+1} = reshape(ParamList,[],1);
                    end
                    obj.S(Param_Idx(1:length(Param)),1) = Param;
                    obj.S_scale(Param_Idx(1:length(Param)),1) = Param;
                    obj.S_ub(Param_Idx,1) = Model('Parameter_UpperBound');
                    obj.S_lb(Param_Idx,1) = Model('Parameter_LowerBound');
                    obj.s_m(Param_Idx,1) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx(1:length(Param)),1) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param,Param_Idx(1:length(Param)),'UniformOutput',false);
                else
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(ParamDescription));
                    for ii = 1:length(bShared)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(cellfun(@(S)strcmp(S,ParamDescription{ii}),obj.s_ParamDescription),1);
                            if isempty(s_Idx)
                                obj.s = obj.s+1; % add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; % add the parameter at the end
                            obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    
                    BranchList = Branch;
                    PrimaryBranch = abs(Branch(1));
                    obj.b_FittingDescription{PrimaryBranch}{end+1} = Model('Model_Description');
                    obj.b_Pdrop{PrimaryBranch}{end+1} = @(q,s)Model({'Pdrop','dPdQ','dPdS'},q,s);
                    obj.b_dPdQ{PrimaryBranch}{end+1} = @(q,s)Model('dPdQ',q,s);
                    obj.b_dPdS{PrimaryBranch}{end+1} = @(q,s)Model('dPdS',q,s);
                    obj.b_Qidx{PrimaryBranch}{end+1} = reshape(BranchList,[],1);
                    obj.b_Sidx{PrimaryBranch}{end+1} = reshape(Param_Idx,[],1);
                    
                    obj.S(Param_Idx(1:length(Param)),1) = Param;
                    obj.S_ub(Param_Idx,1) = Model('Parameter_UpperBound');
                    obj.S_lb(Param_Idx,1) = Model('Parameter_LowerBound');
                    obj.s_m(Param_Idx,1) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx(1:length(Param)),1) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param, Param_Idx(1:length(Param)),'UniformOutput',false);
                end
            elseif ischar(Model)
                Model = DuctNetwork.FittingDatabase(Model);
                obj.AddFitting(Branch,Model,Param);
            elseif iscell(Model)
                cellfun(@(a,b)obj.AddFitting(Branch,a,b),Model,Param);
            end
        end
        
        function Param_Idx = AddParameter(obj, Description, ParamValue, ParamLowerBound, ParamUpperBound, varargin)
            if iscell(Description)
                VAR=cellfun(@(a)num2cell(a),varargin,'UniformOutput',false);
                Param_Idx = cellfun(@(a,b,c,d,varargin)obj.AddParameter(a,b,c,d,varargin{:}),Description, num2cell(ParamValue),num2cell(ParamLowerBound),num2cell(ParamUpperBound),VAR{:});
            elseif isa(Description,'function_handle')
                Model = Description;
                Description = Model('Parameter_Description');
                Is_Identified_Parameter = Model('Is_Identified_Parameter');
                Param_Idx = AddParameter(obj, Description, ParamValue, ParamLowerBound, ParamUpperBound, Is_Identified_Parameter);
            elseif ischar(Description)
                obj.s = obj.s + 1;
                Param_Idx = obj.s;
                obj.s_ParamDescription{Param_Idx,1} = Description;
                obj.S(Param_Idx,1) = ParamValue;
                obj.S_scale(Param_Idx,1) = ParamValue;
                obj.S_ub(Param_Idx,1) = ParamUpperBound;
                obj.S_lb(Param_Idx,1) = ParamLowerBound;
                if nargin>=6
                    obj.s_m(Param_Idx,1) = varargin{1};
                    obj.s_MultiS{Param_Idx,1} = ParamValue*ones(1,varargin{1});
                else
                    obj.s_m(Param_Idx,1) = 0;
                    obj.s_MultiS{Param_Idx,1} = [];
                end
            end
        end
        
        function varargout = SetFlowLimit(obj, q_lb, q_ub)
            x1 = obj.U'*q_lb;
            x2 = obj.U'*q_ub;
            obj.X_ub = x1.*(x1>x2)+x2.*(x2>x1);
            obj.X_lb = x1.*(x1<x2)+x2.*(x2<x1);
            if nargout>0
                varargout{1} = [obj.X_lb;obj.X_ub];
            end
        end
        
        function varargout = BranchPressureDrop(obj, Branch_Idx, X, S) % [dP, dPdX, dPdS]
            N = length(obj.b_FittingDescription{Branch_Idx});
            q = obj.U*X;
            dP = zeros(N,1);dPdQ = zeros(N,obj.b);dPdX = zeros(N,obj.t);dPdS = zeros(N,obj.s);
            if nargout>=3
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    [dP(k), dPdQ(k,idx_Q), dPdS(k,idx_S)] = obj.b_Pdrop{Branch_Idx}{k}(dir.*q(idx_Q),S(idx_S));
                    dP(k) = dir(1)*dP(k);
                    dPdX(k,:) = (dir(1)*dPdQ(k,idx_Q).*dir')*obj.U(idx_Q,:);
                    dPdS(k,:) = dir(1)*dPdS(k,:);
                end
                varargout{1} = sum(dP,1);
                varargout{2} = sum(dPdX,1);
                varargout{3} = sum(dPdS,1);
            elseif nargout==2
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    [dP(k), dPdQ(k,idx_Q)] = obj.b_Pdrop{Branch_Idx}{k}(dir.*q(idx_Q),S(idx_S));
                    %                 Jac2 = [dPdQ(k,idx_Q),dPdS(k,idx_S)];
                    %                 Jac1 = jacobianest(@(x)obj.b_Pdrop{Branch_Idx}{k}(x(1:length(idx_Q)),x(length(idx_Q)+1:end)),[Q(idx_Q);S(idx_S)]);
                    %                 if any(abs(Jac2-Jac1)./Jac1>1e-8)
                    %                     disp(['Jacobian Estimation Error in',obj.b_FittingDescription{Branch_Idx}{k}]);
                    %                     disp([find(abs(Jac2-Jac1)./Jac1>1e-6), max(abs(Jac2-Jac1)./Jac1), dir]);
                    %                 end
                    dP(k) = dir(1)*dP(k);
                    dPdX(k,:) = (dir(1)*dPdQ(k,idx_Q).*dir')*obj.U(idx_Q,:);
                end
                varargout{1} = sum(dP,1);
                varargout{2} = sum(dPdX,1);
            elseif nargout==1
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    dP(k) = obj.b_Pdrop{Branch_Idx}{k}(dir.*q(idx_Q),S(idx_S));
                    dP(k) = dir(1)*dP(k);
                end
                varargout{1} = sum(dP,1);
            end
        end
        
        function [Q,P] = X2QP(obj,X,S_value)
            if nargin<3
                S_value = obj.S;
            end
            dP=cell2mat(cellfun(@(X_i)arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,X_i,S_value),(1:obj.b)'),num2cell(X,1),'UniformOutput',false));
            Q = obj.U*X;
            P = (obj.A')\dP;
        end
        
        function varargout = res_Sim(obj, X, S)
            % Input: X is obj.t dimensional vector of internal state
            %        S is obj.s dimensional vector of system parameter used in simulation
            % Output: e is the residual of State Equations of size obj.t by 1
            %         dedX is the Jacobian of e wrt internal state X of size obj.t by obj.t
            %         dedS is the Jacobian of e wrt system parameter S of size obj.t by obj.s
            X = real(X);
            if nargin==2
                S = obj.S;
            end
            if nargout>=3
                [dP, dPdX, dPdS]=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP); dPdX = cell2mat(dPdX); dPdS = cell2mat(dPdS);
                e = obj.U'*dP*obj.StateResidualControlScale;
                dedX = obj.U'*dPdX*obj.StateResidualControlScale;
                dedS = obj.U'*dPdS*obj.StateResidualControlScale;
                varargout = cell(1,nargout);
                varargout{1}=e; varargout{2}=dedX; varargout{3}=dedS;
                varargout{4}=dP;varargout{5} = dPdX;varargout{6} = dPdS;
                if nargout>=7
                    dXdS = -dedX\dedS;varargout{7} = dXdS;
                end
            elseif nargout==2
                [dP, dPdX]=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP); dPdX = cell2mat(dPdX);
                e = obj.U'*dP*obj.StateResidualControlScale;
                dedX = obj.U'*dPdX*obj.StateResidualControlScale;
                varargout={e,dedX};
            elseif nargout==1
                dP=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP);
                e = obj.U'*dP*obj.StateResidualControlScale;
                varargout={e};
            end
        end
        
        function [X, Q, P] = Sim(obj, Param_Value, Param_Idx)
            if nargin==1
                S_value = obj.S;
            elseif nargin==2
                S_value = Param_Value;
            elseif nargin>=3
                S_value = obj.S;
                if iscell(Param_Idx), Param_Idx = cellfun(@(Str) find(strcmp(Str,obj.s_ParamDescription)),Param_Idx); end
                S_value(Param_Idx) = Param_Value;
            end
            if verLessThan('matlab','9')
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-6,'TolX',1e-12,...
                    'MaxIter',obj.t*20,...
                    'Jacobian','on','DerivativeCheck','off',...
                    'FinDiffType','central','FinDiffRelStep',1e-10,...
                    'OutputFcn',[]);
            else
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'FunctionTolerance',1e-6,'StepTolerance',1e-12,...
                    'MaxIterations',obj.t*20,'OutputFcn',[],...
                    'SpecifyObjectiveGradient',true,'CheckGradients',false,...
                    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-10);
            end
            
            exitflag = -1;resnorm=10;
            X0 = randn(obj.t,1);
            k =obj.Max_Trails;
            resnorm_best = inf;
            while exitflag<=0 || resnorm>1e-6
                k = k - 1;
                obj.n_trail = obj.n_trail+1;
                [X,resnorm,~,exitflag] = lsqnonlin(@(x) obj.res_Sim(x,S_value),X0,obj.X_lb,obj.X_ub,options);
                if resnorm_best>resnorm
                    resnorm_best = resnorm;
                    X_best = X;
                end
                if k>0
                    X = X_best;
                    break;
                end
                X0 = randn(obj.t,1);
            end
            obj.X = X;
            [Q,P] = X2QP(obj,X,S_value);
            obj.Q = Q;
            obj.P = P;
        end
        
        function [r, drdVar] = res_SimCon(obj, Var, num_s, X0, V, D, P0, S0, UndeterminedParam_Idx)
            % X = [internal_state2; UndeterminedParam]
            % r = [res_Sim, res_PressureConstrain]
            if num_s>0
                X_ = X0 + V*Var(1:num_s);
            else
                X_ = X0;
            end
            S_ = S0;
            S_(UndeterminedParam_Idx) = Var(num_s+1:end);
            [e, dedX, dedS, dP, dPdX, dPdS] = obj.res_Sim(X_, S_);
            if isempty(P0) %% no pressure constrain
                r = e;
                drdVar = [dedX*V,dedS(:,UndeterminedParam_Idx)];
            else
                eP = D*dP-P0;
                r = [e; eP];
                deds = dedX*V;
                dedUndetParam = dedS(:,UndeterminedParam_Idx);
                dePds = D*dPdX*V;
                dePdUndetParam = D*dPdS(:,UndeterminedParam_Idx);
                drdVar = [deds,dedUndetParam; dePds, dePdUndetParam];
            end
        end
        
        function [S_val, X, Q, P] = SimCon(obj, Param_Value, Param_Idx, UndeterminedParam_Idx, Q_Constrain_val, Q_Constrain_Idx, P_Constrain_val, P_Constrain_Idx)
            % use one solver and minimum unknown variable to solve sim constrained problem
            % undetermined variable = [remained_state; UndeterminedParam]
            % residual = [res_Sim, res_PressureConstrain]
            
            C = sparse(1:length(Q_Constrain_Idx),Q_Constrain_Idx,1,length(Q_Constrain_Idx),obj.b);
            X0 = (C*obj.U)\reshape(Q_Constrain_val,[],1);
            V = null(C*obj.U);
            num_s = obj.t - rank(C*obj.U);
            if isempty(P_Constrain_Idx)
                D = [];
                P0 = [];
            else
                D = sparse(1:length(P_Constrain_Idx),P_Constrain_Idx,1,length(P_Constrain_Idx),obj.n)*((obj.A')\eye(obj.b));
                P0 = reshape(P_Constrain_val,[],1);
            end
            S0 = obj.S;
            S0(Param_Idx) = Param_Value;
            
            if verLessThan('matlab','9')
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-12,'TolX',1e-12,...
                    'MaxIter',5*(num_s+length(UndeterminedParam_Idx)),...
                    'Jacobian','on',...
                    'OutputFcn',[]);
            else
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12,...
                    'OptimalityTolerance',1e-12,...
                    'MaxIterations',5*(num_s+length(UndeterminedParam_Idx)),'OutputFcn',[],...
                    'SpecifyObjectiveGradient',true);
            end
            Var0 = [rand(num_s,1);S0(UndeterminedParam_Idx)];
            Var_ub = [ 1e12*ones(num_s,1);obj.S_ub(UndeterminedParam_Idx)];
            Var_lb = [-1e12*ones(num_s,1);obj.S_lb(UndeterminedParam_Idx)];
            f = @(x)obj.res_SimCon(x, num_s, X0, V, D, P0, S0, UndeterminedParam_Idx);
            [Var, resnorm] = lsqnonlin(f,Var0,Var_lb,Var_ub,options);
            X = X0 + V*Var(1:num_s);
            S_val = Var(1+num_s:end);
            S_ = S0; S_(UndeterminedParam_Idx) = S_val;
            if resnorm>1
                disp('fail to satisfy the constrain');
                [X,Q,P] = obj.Sim(S_);
            else
                [X,Q,P] = obj.Sim(S_);
                %[Q,P] = obj.X2QP(X,S_);
            end
        end
        
        function [R, dRdVar] = res_MultiSimCon(obj, Var, cnfg_VarIdx, cnfg_num_s, cnfg_X0, cnfg_V, cnfg_D, cnfg_P0, cnfg_S0, cnfg_UndeterminedParam_Idx)
            cnfg_x = cellfun(@(idx)Var(idx),cnfg_VarIdx,'UniformOutput',false);
            [cnfg_r, cnfg_drdVar]=cellfun(@obj.res_SimCon,cnfg_x,cnfg_num_s, cnfg_X0, cnfg_V, cnfg_D, cnfg_P0, cnfg_S0, cnfg_UndeterminedParam_Idx,'UniformOutput',false);
            R = cell2mat(cnfg_r);
            cnfg_VarMap = cellfun(@(idx)sparse(1:length(idx),idx,1,length(idx),length(Var)),cnfg_VarIdx,'UniformOutput',false);
            cnfg_dRdVar = cellfun(@mtimes,cnfg_drdVar,cnfg_VarMap,'UniformOutput',false);
            dRdVar = cell2mat(cnfg_dRdVar);
        end
        
        function [cnfg_UndetParam_Val, cnfg_X, cnfg_Q, cnfg_P] = MultiSimCon(obj, cnfg_Param_Value, cnfg_Param_Idx, UndetParam_Idx, cnfg_UndetParam_Ord, cnfg_Q_Constrain_val, cnfg_Q_Constrain_Idx, cnfg_P_Constrain_val, cnfg_P_Constrain_Idx)
            % Prepare cnfg_VarIdx, cnfg_num_s, cnfg_X0, cnfg_V, cnfg_D, cnfg_P0, cnfg_S0, cnfg_UndeterminedParam_Idx
            cnfg = length(cnfg_Q_Constrain_Idx);
            num_UndetParam_Idx = length(UndetParam_Idx);
            num_UndetParam_Ord = max(cat(2,cnfg_UndetParam_Ord{:}),[],2);
            num_UndetParam = sum(num_UndetParam_Ord);
            UndetParam_Loc = [0;cumsum(num_UndetParam_Ord(1:end-1))];
            
            cnfg_C = cellfun(@(idx)sparse(1:length(idx),idx,1,length(idx),obj.b),cnfg_Q_Constrain_Idx,'UniformOutput',false);
            num_s = cellfun(@(C)obj.t - rank(C*obj.U),cnfg_C,'UniformOutput',true);
            sum_num_s = sum(num_s);
            cnfg_VarIdx_1 = mat2cell((1:sum_num_s)',num_s,1);
            cnfg_VarIdx_2 = cellfun(@(order)sum_num_s+UndetParam_Loc+order,cnfg_UndetParam_Ord,'UniformOutput',false);
            cnfg_VarIdx = cellfun(@(a,b)[a;b],cnfg_VarIdx_1,cnfg_VarIdx_2,'UniformOutput',false);%cnfg_VarIdx
            
            cnfg_num_s = num2cell(num_s);
            cnfg_X0 = cellfun(@(C,Q_Constrain_val)(C*obj.U)\reshape(Q_Constrain_val,[],1),cnfg_C,cnfg_Q_Constrain_val,'UniformOutput',false);
            cnfg_V = cellfun(@(C)null(C*obj.U),cnfg_C,'UniformOutput',false);
            if isempty(cnfg_P_Constrain_Idx)
                cnfg_D = repmat({[]},cnfg,1);
                cnfg_P0 = repmat({[]},cnfg,1);
            else
                cnfg_D = cellfun(@(idx)sparse(1:length(idx),idx,1,length(idx),obj.n)*((obj.A')\eye(obj.b)),cnfg_P_Constrain_Idx,'UniformOutput',false);
                cnfg_P0 = cellfun(@(val)reshape(val,[],1),cnfg_P_Constrain_val,'UniformOutput',false);
            end
            if iscell(cnfg_Param_Value)
                cnfg_S0 = repmat({reshape(obj.S,[],1)},cnfg,1);
                for kk = length(cnfg_Param_Value)
                    if isempty(cnfg_Param_Idx{kk})
                        cnfg_S0{kk}=cnfg_Param_Value{kk};
                    else
                        cnfg_S0{kk}(cnfg_Param_Idx{kk})=cnfg_Param_Value{kk};
                    end
                end
            elseif isempty(cnfg_Param_Value)
                cnfg_S0 = repmat({reshape(obj.S,[],1)},cnfg,1);
            else %numeric array
                S0 = obj.S;
                S0(cnfg_Param_Idx)=cnfg_Param_Value;
                cnfg_S0 = repmat({reshape(S0,[],1)},cnfg,1);
            end
            cnfg_UndetParam_Idx = repmat({UndetParam_Idx},cnfg,1);
            % Initialization
            Var0 = rand(sum_num_s + num_UndetParam,1);
            f = @(x)obj.res_MultiSimCon(x, cnfg_VarIdx, cnfg_num_s, cnfg_X0, cnfg_V, cnfg_D, cnfg_P0, cnfg_S0, cnfg_UndetParam_Idx);
            if length(f(Var0))<length(Var0)
                Var_ub = [];
                Var_lb = [];
            else
                Var_ub = [ 1e12*ones(sum_num_s,1);repelem(obj.S_ub(UndetParam_Idx),num_UndetParam_Ord);];
                Var_lb = [-1e12*ones(sum_num_s,1);repelem(obj.S_lb(UndetParam_Idx),num_UndetParam_Ord);];
            end
            if verLessThan('matlab','9')
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-12,'TolX',1e-12,...
                    'MaxIter',15*length(Var0),...
                    'Jacobian','on',...
                    'OutputFcn',[]);
            else
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12,...
                    'OptimalityTolerance',1e-12,...
                    'MaxIterations',15*length(Var0),'OutputFcn',[],...
                    'SpecifyObjectiveGradient',true);
            end
            % Solve for undetermined variables by optimization
            [Var, resnorm] = lsqnonlin(f,Var0,Var_lb,Var_ub,options);
            
            % Analyze the result
            % cnfg_s = mat2cell(Var(1:sum_num_s),num_s,1);
            % cnfg_X = cellfun(@(X0,V,s)X0 + V*s,cnfg_X0,cnfg_V,cnfg_s,'UniformOutput',false);
            UndetParam_Val = Var((1:num_UndetParam)+sum_num_s);
            cnfg_UndetParam_Val = cellfun(@(order)UndetParam_Val(UndetParam_Loc+order),cnfg_UndetParam_Ord,'UniformOutput',false);
            cnfg_S = cnfg_S0;
            for kk = 1:cnfg
                cnfg_S{kk}(UndetParam_Idx) = cnfg_UndetParam_Val{kk};
            end
            DispState = obj.Display;
            obj.Display = 'none';
            [cnfg_X,cnfg_Q,cnfg_P] = cellfun(@obj.Sim,cnfg_S,'UniformOutput',false);
            obj.Display = DispState;
        end
        
        function [X, S_val, Q, P] = SimConstrained(obj, Param_Value, Param_Idx, UndeterminedParam_Idx, Q_Constrain_val, Q_Constrain_Idx, P_Constrain_val, P_Constrain_Idx)
            % Param_Value is array to specify parameters
            % Param_Idx is array or cell, which indicates the index of parameters in obj.S that uses Param_Value, others uses obj.S
            % UndeterminedParam_Idx is array or cell, which indicates the index of parameters that need to be determined
            
            % Q is cell of size b by 1, scalar value specifies the target value, set as [] if no specification on it
            % P is cell of size n by 1, scalar value specifies the target value, set as [] if no specification on it
            if iscell(Param_Idx)
                Param_Idx = cellfun(@(Str) find(strcmp(Str,obj.s_ParamDescription)),Param_Idx);
            end
            if iscell(UndeterminedParam_Idx)
                UndeterminedParam_Idx = cellfun(@(Str) find(strcmp(Str,obj.s_ParamDescription)),UndeterminedParam_Idx);
            end
            Backup = obj.BackupParameterIdentification();
            cnfg_s_mIdx = zeros(1,obj.s);
            cnfg_s_mIdx(UndeterminedParam_Idx) = 1;
            obj.S(Param_Idx) = Param_Value;
            
            proc_ob_SensorLocation = [Q_Constrain_Idx(:);P_Constrain_Idx(:)]';
            ob_SensorType = [false(size(Q_Constrain_Idx(:)));true(size(P_Constrain_Idx(:)))];
            proc_ob_Measurements = [Q_Constrain_val;P_Constrain_val]';
            Q_Uncertainty = 1e-3;
            P_Uncertainty = 2;
            obj.ob_Uncertainty = (~ob_SensorType)*Q_Uncertainty + ob_SensorType*P_Uncertainty;
            [MultiS, cnfg_X] = obj.Identification(cnfg_s_mIdx, 1, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, obj.ob_Uncertainty);
            X = cnfg_X{1};
            S_value = obj.S;
            S_value(sort(UndeterminedParam_Idx)) = MultiS; % MultiS are always in ascending order from the first parameter to last
            [Q,P] = X2QP(obj,X,S_value);
            S_val = S_value(UndeterminedParam_Idx);
            obj.RecoverParameterIdentification(Backup);
        end
        
        function [IdentifiedVal, X_final, Q_final, P_final] = SimConstrained_CC(obj, ModifiedParam_Idx, ModifiedParamParam_Value, Constrain_QIdx, ConstrainVal_Q, Constrain_PIdx, ConstrainVal_P, UndetParamIdx, SensorAccuracy)
            % Constrain_QIdx/Constrain_PIdx is index in Q_list/P_list
            
            QIdx = Constrain_QIdx;
            MeasureData_Q = ConstrainVal_Q;
            PIdx = Constrain_PIdx;
            MeasureData_P = ConstrainVal_P;
            
            f = @(x)Res_SimConstrained(x, obj, QIdx, MeasureData_Q, PIdx, MeasureData_P, UndetParamIdx, SensorAccuracy);
            x0 = rand(length(UndetParamIdx),1);
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,...
                'CheckGradients',false,'FiniteDifferenceType','central',...
                'Display','off',...
                'FunctionTolerance',1e-12,...
                'StepTolerance',1e-8,...
                'MaxIterations',400);
            IdentifiedVal = lsqnonlin(f,x0,zeros(size(x0)),[],options);
            
            obj.S(ModifiedParam_Idx) = ModifiedParamParam_Value;
            S = obj.S;
            S(UndetParamIdx) = IdentifiedVal;
            [X_final, Q_final, P_final] = obj.Sim(S);
            
            function [residual, Jacobian] = Res_SimConstrained(X, obj, QIdx, MeasureData_Q, PIdx, MeasureData_P, UndetParamIdx, SensorAccuracy)
                S_P = obj.S;
                S_P(UndetParamIdx) = X;
                
                [X_internal, Q, P] = obj.Sim(S_P);
                [~,~,~,~,ddPdX,ddPdS,dXdS] = obj.res_Sim(X_internal, S_P);
                dQdS = obj.U*dXdS;
                dPdS = (obj.A*obj.A')\obj.A*(ddPdX*dXdS+ddPdS);
                
                PredictData_Q = Q(QIdx);
                PredictData_P = P(PIdx);
                Jacobian_Q = dQdS(QIdx,UndetParamIdx);
                Jacobian_P = dPdS(PIdx,UndetParamIdx);
                
                residual_Q = (PredictData_Q - MeasureData_Q)/SensorAccuracy(1);
                residual_P = (PredictData_P - MeasureData_P)/SensorAccuracy(2);
                residual = [residual_Q;residual_P];
                Jacobian = [Jacobian_Q/SensorAccuracy(1);Jacobian_P/SensorAccuracy(2)];
            end
        end
        
        function Q = SetDamperAndReadFlow(obj, theta, DamperID)
            if ~isempty(theta)
                obj.S(obj.d_Sidx(DamperID))=theta;
            end
            if nargout>0
                obj.Sim();
                Q = reshape(obj.Q(obj.gdr_Bidx(DamperID)),1,[]);
            end
        end
        
        function proc_ob_Measurements = Measurements(obj, ctrl_ParamIdx, proc_ctrl_ParamValue, ob_SensorType, ob_Uncertainty, proc_ob_SensorPosition)
            obj.proc = size(proc_ctrl_ParamValue,1);
            obj.ob = length(ob_SensorType);
            if isreal(ob_Uncertainty)
                ob_UncertaintyArray = ob_Uncertainty;
                ob_Uncertainty = cell(size(ob_Uncertainty));
                for ii=1:obj.ob
                    ob_Uncertainty{ii} = @(x) x+randn()*ob_UncertaintyArray(ii);
                end
            end
            proc_ob_Measurements  = zeros(obj.proc,obj.ob);
            for ii = 1:obj.proc
                obj.Sim(proc_ctrl_ParamValue(ii,:),ctrl_ParamIdx);
                for jj = 1:obj.ob
                    switch ob_SensorType{jj}
                        case {'P','Pressure'}
                            proc_ob_Measurements(ii,jj) = ob_Uncertainty{jj}(obj.P(proc_ob_SensorPosition(ii,jj)));
                        case {'q','Flowrate'}
                            proc_ob_Measurements(ii,jj) = ob_Uncertainty{jj}(obj.Q(proc_ob_SensorPosition(ii,jj)));
                    end
                end
            end
        end
        
        function proc_ob_Measurements = Chen_Method_Measurement(obj, Damper_sIdx, DamperPositions, Flow_bIdx, Pressure_nIdx, Flow_Uncertainty, Pressure_Uncertainty)
            % Damper_sIdx vector indicating damper parameters' index in S, length(Damper_sIdx) = d
            % DamperPositions vector indicating damper positions, length(DamperPositions) = m
            % Flow_bIdx vector indicating flow meter locations in branch ID, length(Flow_bIdx) = d
            % Pressure_nIdx vector indicating pressure sensor locations in node ID, length(Pressure_nIdx) = d
            % Flow_Uncertainty is the uncertainty superposed on flow meter;
            % Pressure_Uncertainty is the uncertainty superposed on pressure sensor'
            
            m = length(DamperPositions);
            d = length(Damper_sIdx);
            ctrl_ParamIdx = Damper_sIdx(:);
            proc_ctrl_ParamValue = full(sparse(1:m*d,reshape(ones(m,1)*(1:d),1,[]),reshape(repmat(DamperPositions(:)-DamperPositions(1),1,d),1,[]),m*d,d)) + DamperPositions(1);
            ob_SensorType = {'q';'P'};
            ob_Uncertainty = [Flow_Uncertainty;Pressure_Uncertainty];
            proc_ob_SensorPosition = [reshape(ones(m,1)*Flow_bIdx(:)',[],1),reshape(ones(m,1)*Pressure_nIdx(:)',[],1)];
            proc_ob_Measurements = obj.Measurements(ctrl_ParamIdx, proc_ctrl_ParamValue, ob_SensorType, ob_Uncertainty, proc_ob_SensorPosition);
        end
        
        function varargout = res_Identification(obj, Tot_r)
            % Input: Tot_r is the all variables in identification, Tot_r = [cnfg_X(:);[s_MultiS{:}]']
            %        cnfg_X is the internal states of each configuration
            %        number matrix of size t by cnfg s_MultiS is all identified parameters packed in cell array of size s by 1 in which length of number array in each cell is s_m
            % Output: Tot_e is the residual of all identification equations, Tot_e = [cnfg_e{:};proc_e{:}]
            %         cnfg_e is the residual of state equations for each configuration
            %         number array of size t by cnfg (folded horizontally) proc_e is the residual of measurements in each procedure
            %         number array of size ob by proc (folded horizontally) dTot_edTot_r
            %             disp(reshape(Tot_r,1,[]))
            cnfg_X = num2cell(reshape(Tot_r(1:obj.cnfg*obj.t),obj.t,obj.cnfg),1)';
            cnfg_S = repmat({obj.S},obj.cnfg,1);
            IdentSIdx = obj.s_m>0;
            for ii = 1:obj.cnfg
                cnfg_S{ii}(IdentSIdx)=Tot_r(obj.cnfg*obj.t+obj.cnfg_MultiSIdx{ii}(IdentSIdx));
                %Each element of cnfg_MultiSIdx contains the idexing of identified parameters in the MultiS vector
            end
            
            [cnfg_e,cnfg_dedX,cnfg_dedS, cnfg_dP, cnfg_dPdX, cnfg_dPdS] = cellfun(@obj.res_Sim,cnfg_X,cnfg_S,'UniformOutput',false);
            cnfg_QP = cellfun(@(X,dP)[obj.U*X; (obj.A')\dP],cnfg_X,cnfg_dP,'UniformOutput',false);
            proc_mu = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_QP(obj.proc_ConfigIdx,:),'UniformOutput',false);
            proc_e = cellfun(@(mu,Z)rdivide(mu-Z,obj.ob_Uncertainty),proc_mu,num2cell(obj.proc_Measurements',1)','UniformOutput',false);
            
            Tot_e = cell2mat([cnfg_e;proc_e]);
            varargout{1} = Tot_e;
            if nargout>1
                cnfg_dQPdX = cellfun(@(dPdX)[obj.U;(obj.A')\dPdX],cnfg_dPdX,'UniformOutput',false);
                cnfg_dQPdS = cellfun(@(dPdS)[zeros(obj.b,obj.s);(obj.A')\dPdS],cnfg_dPdS,'UniformOutput',false);
                %cnfg_dQPdS = cellfun(@(dPdS,dPdX,dXdS)[obj.U*dXdS;(obj.A')\(dPdS+dPdX*dXdS)],cnfg_dPdS,cnfg_dPdX,cnfg_dXdS,'UniformOutput',false);%
                proc_dmudX = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_dQPdX(obj.proc_ConfigIdx,:),'UniformOutput',false);
                proc_dmudS = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_dQPdS(obj.proc_ConfigIdx,:),'UniformOutput',false);
                proc_dedX = cellfun(@(dmudX)rdivide(dmudX,repmat(obj.ob_Uncertainty,1,obj.t)),proc_dmudX,'UniformOutput',false);
                proc_dedS = cellfun(@(dmudS)rdivide(dmudS,repmat(obj.ob_Uncertainty,1,obj.s)),proc_dmudS,'UniformOutput',false);
                %                 cnfg_dSdMultiS = cellfun(@(MultiSIdx)sparse(IdentSIdx,MultiSIdx,1,obj.s,sum(obj.s_m)),obj.cnfg_MultiSIdx,'UniformOutput',false);
                %                 cnfg_dedMultiS = cellfun(@mtimes,cnfg_dedS,cnfg_dSdMultiS,'UniformOutput',false);
                dcnfg_edcnfg_X = blkdiag(cnfg_dedX{:});
                dcnfg_edMultiS = zeros(obj.t*obj.cnfg,sum(obj.s_m));
                for ii = 1:obj.cnfg
                    dcnfg_edMultiS(obj.t*(ii-1)+(1:obj.t),obj.cnfg_MultiSIdx{ii}(IdentSIdx))=cnfg_dedS{ii}(:,IdentSIdx);
                end
                dproc_edcnfg_X = zeros(obj.ob*obj.proc,obj.t*obj.cnfg);
                dproc_edMultiS = zeros(obj.ob*obj.proc,sum(obj.s_m));
                for ii = 1:obj.proc
                    RowIdx = obj.ob*(ii-1)+(1:obj.ob);
                    dproc_edcnfg_X(RowIdx,obj.t*(obj.proc_ConfigIdx(ii)-1)+(1:obj.t)) = proc_dedX{ii};
                    dproc_edMultiS(RowIdx,obj.cnfg_MultiSIdx{obj.proc_ConfigIdx(ii)}(IdentSIdx)) = proc_dedS{ii}(:,IdentSIdx);
                end
                dTot_edTot_r = [dcnfg_edcnfg_X,dcnfg_edMultiS;dproc_edcnfg_X,dproc_edMultiS];
                varargout{2} = dTot_edTot_r;
            end
        end
        
        function varargout = res_Identification_Cascaded(obj, Tot_r)
            % Input: Tot_r is the all variables in identification, Tot_r = [s_MultiS{:}]'
            %        s_MultiS is all identified parameters packed in cell array of size s by 1 in which length of number array in each cell is s_m
            % Output: Tot_e is the residual of all identification equations, Tot_e = [proc_e{:}]
            %        number array of size ob by proc (folded horizontally) proc_e is the residual of measurements in each procedure
            %         dTot_edTot_r is the jacobian of the function with size ob*proc by sum(s_m)
            
            %        cnfg_X is the internal states of each configuration
            %        cnfg_e is the residual of state equations for each configuration
            
            cnfg_S = repmat({obj.S},obj.cnfg,1);
            IdentSIdx = obj.s_m>0;
            for ii = 1:obj.cnfg
                cnfg_S{ii}(IdentSIdx)=Tot_r(obj.cnfg_MultiSIdx{ii}(IdentSIdx));
                %Each element of cnfg_MultiSIdx contains the idexing of identified parameters in the MultiS vector
            end
            backup.Display = obj.Display;
            obj.Display = 'none';
            [cnfg_X, cnfg_Q, cnfg_P] = cellfun(@obj.Sim,cnfg_S,'UniformOutput',false);
            obj.Display = backup.Display;
            [~,~,~,~, cnfg_dPdX, cnfg_dPdS, cnfg_dXdS] = cellfun(@obj.res_Sim,cnfg_X,cnfg_S,'UniformOutput',false);
            %             cnfg_QP = cellfun(@(X,dP)[obj.U*X; (obj.A')\dP],cnfg_X,cnfg_dP,'UniformOutput',false);
            cnfg_QP = cellfun(@vertcat,cnfg_Q,cnfg_P,'UniformOutput',false);
            proc_mu = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_QP(obj.proc_ConfigIdx,:),'UniformOutput',false);
            proc_e = cellfun(@(mu,Z)rdivide(mu-Z,obj.ob_Uncertainty),proc_mu,num2cell(obj.proc_Measurements',1)','UniformOutput',false);
            Tot_e = cell2mat(proc_e);
            varargout{1} = Tot_e;
            if nargout>1
                %cnfg_dQPdX = cellfun(@(dPdX)[obj.U;(obj.A')\dPdX],cnfg_dPdX,'UniformOutput',false);
                %cnfg_dQPdS = cellfun(@(dPdS)[zeros(obj.b,obj.s);(obj.A')\dPdS],cnfg_dPdS,'UniformOutput',false);
                cnfg_dQPdS = cellfun(@(dPdS,dPdX,dXdS)[obj.U*dXdS;(obj.A')\(dPdS+dPdX*dXdS)],cnfg_dPdS,cnfg_dPdX,cnfg_dXdS,'UniformOutput',false);%
                %proc_dmudX = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_dQPdX(obj.proc_ConfigIdx,:),'UniformOutput',false);
                proc_dmudS = cellfun(@mtimes,obj.proc_ObPosMatrix,cnfg_dQPdS(obj.proc_ConfigIdx,:),'UniformOutput',false);
                %proc_dedX = cellfun(@(dmudX)rdivide(dmudX,repmat(obj.ob_Uncertainty,1,obj.t)),proc_dmudX,'UniformOutput',false);
                proc_dedS = cellfun(@(dmudS)rdivide(dmudS,repmat(obj.ob_Uncertainty,1,obj.s)),proc_dmudS,'UniformOutput',false);
                %                 cnfg_dSdMultiS = cellfun(@(MultiSIdx)sparse(IdentSIdx,MultiSIdx,1,obj.s,sum(obj.s_m)),obj.cnfg_MultiSIdx,'UniformOutput',false);
                %                 cnfg_dedMultiS = cellfun(@mtimes,cnfg_dedS,cnfg_dSdMultiS,'UniformOutput',false);
                %                 dcnfg_edcnfg_X = blkdiag(cnfg_dedX{:});
                %                 dcnfg_edMultiS = zeros(obj.t*obj.cnfg,sum(obj.s_m));
                %                 for ii = 1:obj.cnfg
                %                     dcnfg_edMultiS(obj.t*(ii-1)+(1:obj.t),obj.cnfg_MultiSIdx{ii}(IdentSIdx))=cnfg_dedS{ii}(:,IdentSIdx);
                %                 end
                %                 dproc_edcnfg_X = zeros(obj.ob*obj.proc,obj.t*obj.cnfg);
                dproc_edMultiS = zeros(obj.ob*obj.proc,sum(obj.s_m));
                for ii = 1:obj.proc
                    RowIdx = obj.ob*(ii-1)+(1:obj.ob);
                    %                     dproc_edcnfg_X(RowIdx,obj.cnfg*(obj.proc_ConfigIdx(ii)-1)+(1:obj.t)) = proc_dedX{ii};
                    dproc_edMultiS(RowIdx,obj.cnfg_MultiSIdx{obj.proc_ConfigIdx(ii)}(IdentSIdx)) = proc_dedS{ii}(:,IdentSIdx);
                end
                dTot_edTot_r = dproc_edMultiS;
                varargout{2} = dTot_edTot_r;
                varargout{3} = cnfg_X;
            end
        end
        
        function [MultiS, cnfg_X, R] = Identification(obj, cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty)
            [MultiS, cnfg_X, R] = Identification_Cascaded(obj, cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty);
            %             [MultiS, cnfg_X, R] = Identification_Hybrid(obj, cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty);
        end
        
        function [MultiS, cnfg_X, R] = Identification_Hybrid(obj, cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty)
            % cnfg_s_mIdx is matrix of size cnfg by s, each of the cell indicates the identified parameter and its index in s_MultiS for each configuration
            % proc_ConfigIdx is vector of size proc by 1, indicating each process uses which configuration.
            % proc_ob_Measurements is matrix of size proc by ob, storing the measured data in each process
            % proc_ob_SensorLocation is matrix of size proc by ob, indicating the position of sensors in each process
            % ob_SensorType is vector or cell array of size ob by 1, indicating the type of each sensor, Pressure sensor: 1|'Pressure'|'p', Flowrate sensor: 0|'Flowrate'|'q'
            % ob_Uncertainty is vector of size ob by 1, indicating the uncertainty for each sensor.
            %             Backup = obj.BackupParameterIdentification();
            obj.Identification_Preprocess(cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty);
            if verLessThan('matlab','9')
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'TolFun',1e-12,'TolX',1e-12,...
                    'MaxIter',10,...
                    'Jacobian','on',...
                    'OutputFcn',[]);
            else
                options = optimoptions(@lsqnonlin,'Display',obj.Display,...
                    'Algorithm','trust-region-reflective',...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12,...
                    'OptimalityTolerance',1e-12,...
                    'MaxIterations',10,'OutputFcn',[],...
                    'SpecifyObjectiveGradient',true);
            end
            exitflag = -1;
            n_tot_S = sum(obj.s_m);
            n_tot_X = obj.cnfg*obj.t;
            if isempty(obj.X)
                S0 = ones(n_tot_S,1);
                X0 = ones(n_tot_X,1);
            else
                S0 = cell2mat(arrayfun(@(x,b)repmat(x,b,1),obj.S,obj.s_m,'UniformOutput',false));
                X0 = repmat(obj.X,obj.cnfg,1);
            end
            ub_S = cell2mat(arrayfun(@(x,b)repmat(x,b,1),obj.S_ub,obj.s_m,'UniformOutput',false));
            lb_S = cell2mat(arrayfun(@(x,b)repmat(x,b,1),obj.S_lb,obj.s_m,'UniformOutput',false));
            ub_X = repmat(obj.X_ub,obj.cnfg,1);
            lb_X = repmat(obj.X_lb,obj.cnfg,1);
            obj.StateResidualControlScale = 10;
            curr = inf; resnorm = 1e10;
            while exitflag<=0
                R = lsqnonlin(@obj.res_Identification,[X0;S0],[lb_X;lb_S],[ub_X;ub_S],optimoptions(options,'MaxIterations',max(20,round(-3*log(resnorm)))));
                prev = curr; curr = norm(obj.res_Identification_Cascaded(R(n_tot_X+1:n_tot_S+n_tot_X)))^2;
                if curr<prev % accept the update by res_Identification
                    S0 = R(n_tot_X+1:n_tot_S+n_tot_X);
                    if strcmp(obj.Display,'iter')
                        disp('Accept Update by res_Identification');
                    end
                else
                    if strcmp(obj.Display,'iter')
                        disp('Reject Update by res_Identification');
                    end
                end
                [S0,resnorm] = lsqnonlin(@obj.res_Identification_Cascaded,S0,lb_S,ub_S,options);
                %                 if strcmp(obj.Display,'iter')
                %                     disp(S0');
                %                     system('pause');
                %                 end
                cnfg_S = repmat({obj.S},obj.cnfg,1);IdentSIdx = obj.s_m>0;
                for ii = 1:obj.cnfg
                    cnfg_S{ii}(IdentSIdx)=S0(obj.cnfg_MultiSIdx{ii}(IdentSIdx));
                end
                backup.Display = obj.Display; obj.Display = 'none';
                cnfg_X = cellfun(@obj.Sim,cnfg_S,'UniformOutput',false); obj.Display = backup.Display;
                X0 = cell2mat(cnfg_X);
                
                if resnorm<1e-4
                    exitflag = 1;
                    continue;
                end
            end
            MultiS = S0;
            obj.s_MultiS=mat2cell(MultiS,obj.s_m,1);
        end
        
        function [Backup] = Identification_Preprocess(obj,cnfg_s_mIdx, proc_ConfigIdx, proc_ob_Measurements, proc_ob_SensorLocation, ob_SensorType, ob_Uncertainty)
            Backup.proc_Measurements = obj.proc_Measurements;
            Backup.proc = obj.proc; Backup.ob = obj.ob;
            Backup.proc_ObPosMatrix = obj.proc_ObPosMatrix;
            Backup.proc_ConfigIdx = obj.proc_ConfigIdx;
            Backup.ob_Uncertainty = obj.ob_Uncertainty;
            Backup.cnfg = obj.cnfg;
            Backup.s_m = obj.s_m;
            Backup.cnfg_MultiSIdx = obj.cnfg_MultiSIdx;
            
            obj.proc_Measurements = proc_ob_Measurements;
            [obj.proc, obj.ob] = size(proc_ob_Measurements);
            if ~islogical(ob_SensorType) && ~isnumeric(ob_SensorType)
                ob_SensorStr = ob_SensorType;
                ob_SensorType = true(1,obj.ob);
                for ii = 1:obj.ob
                    switch ob_SensorStr{ii}
                        case {'P','Pressure'}
                            ob_SensorType(ii) = true;
                        case {'q','Flowrate'}
                            ob_SensorType(ii) = false;
                    end
                end
            end
            obj.proc_ObPosMatrix = cellfun(@(ob_SensLoc)sparse(1:obj.ob,ob_SensLoc+reshape(ob_SensorType,1,obj.ob)*obj.b,1,obj.ob,obj.n+obj.b),num2cell(proc_ob_SensorLocation,2),'UniformOutput',false);
            obj.proc_ConfigIdx = proc_ConfigIdx;
            obj.ob_Uncertainty = ob_Uncertainty;
            obj.cnfg = size(cnfg_s_mIdx,1);
            obj.s_m = max(cnfg_s_mIdx,[],1)';
            obj.cnfg_MultiSIdx = cellfun(@(x)(cumsum(obj.s_m)-obj.s_m + x').*(x'>0),num2cell(cnfg_s_mIdx,2),'UniformOutput',false);
        end
        
        function Backup = BackupParameterIdentification(obj)
            Backup.S = obj.S;
            Backup.s_m = obj.s_m;
            Backup.s_MultiS = obj.s_MultiS;
            Backup.cnfg = obj.cnfg;
            Backup.cnfg_MultiSIdx = obj.cnfg_MultiSIdx;
            Backup.ob = obj.ob;
            Backup.ob_Uncertainty = obj.ob_Uncertainty;
            Backup.proc = obj.proc;
            Backup.proc_ConfigIdx = obj.proc_ConfigIdx;
            Backup.proc_ObPosMatrix = obj.proc_ObPosMatrix;
            Backup.proc_Measurements = obj.proc_Measurements;
        end
        
        function RecoverParameterIdentification(obj, Backup)
            obj.S = Backup.S;
            obj.s_m = Backup.s_m;
            obj.s_MultiS = Backup.s_MultiS;
            obj.cnfg = Backup.cnfg;
            obj.cnfg_MultiSIdx = Backup.cnfg_MultiSIdx;
            obj.ob = Backup.ob;
            obj.ob_Uncertainty = Backup.ob_Uncertainty;
            obj.proc = Backup.proc;
            obj.proc_ConfigIdx = Backup.proc_ConfigIdx;
            obj.proc_ObPosMatrix = Backup.proc_ObPosMatrix;
            obj.proc_Measurements = Backup.proc_Measurements;
        end
        
        function [cnfg_s_mIdx, proc_ConfigIdx, proc_ob_SensorLocation, ob_SensorType] = Identification_Parameter_Generator(obj, Damper_sIdx, m, UndetParam_sIdx, Flow_bIdx, Pressure_nIdx)
            % Damper_sIdx vector indicating damper parameters' index in S, length(Damper_sIdx) = d
            % m indicating the number of damper position
            % UndetParam_sIdx vector indicating undetermined parameters' index in S
            % Flow_bIdx vector indicating flow meter locations in branch ID, length(Flow_bIdx) = d
            % Pressure_nIdx vector indicating pressure sensor locations in node ID, length(Pressure_nIdx) = d
            
            d = length(Damper_sIdx);
            cnfg_s_mIdx = repmat(full(sparse(1,[UndetParam_sIdx(:);Damper_sIdx(:)],1,1,obj.s)),1+(m-1)*d,1);
            for ii = 1:d
                cnfg_s_mIdx(1+(m-1)*(ii-1)+(1:m-1),Damper_sIdx(ii)) = (2:m)';
            end
            proc_ConfigIdx = reshape([ones(1,d);reshape(2:(m-1)*d+1,m-1,d)],[],1);
            proc_DamperIdx = reshape(repmat(1:d,m,1),[],1);
            Flow_bIdx = Flow_bIdx(:);
            Pressure_nIdx = Pressure_nIdx(:);
            proc_ob_SensorLocation = [Flow_bIdx(proc_DamperIdx(:)),Pressure_nIdx(proc_DamperIdx(:))];
            ob_SensorType = {'q';'P'};
        end
    end
    
    methods (Static = true)
        function ModelFcn = FittingDatabase(ModelStr)
            ModelFcn=str2func(['DuctNetwork.',ModelStr]);
        end
        
        function varargout = FanQuadratic(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = -s(1);
                        else
                            varargout{ii} = s(1)*(q^2/s(2)^2-1);
                        end
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = 0;
                        else
                            varargout{ii} = s(1)*2*q/s(2)^2;
                        end
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = [-1,0];
                        else
                            varargout{ii} = [q^2/s(2)^2-1, -s(1)*2*q^2/s(2)^3];
                        end
                    case 'Model_Description'
                        varargout{ii}='Fan Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.FanQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Max Pressure(Pa)','Max Flow(m^3/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-3,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = FanQuadratic_kPa(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = -s(1)*1000;
                        else
                            varargout{ii} = s(1)*1000*(q^2/s(2)^2-1);
                        end
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = 0;
                        else
                            varargout{ii} = s(1)*1000*2*q/s(2)^2;
                        end
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = [-1,0]*1000;
                        else
                            varargout{ii} = [q^2/s(2)^2-1, -s(1)*2*q^2/s(2)^3]*1000;
                        end
                    case 'Model_Description'
                        varargout{ii}='Fan Using Quadratic Q-dP Relationship (kPa unit)';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.FanQuadratic_kPa};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Max Pressure(kPa)','Max Flow(m^3/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e3,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = Fan_NVT250(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        [dP, dPdQ, dPdS] = DuctNetwork.Fan(q,s,DuctNetwork.Table_NVT250.Model_Coeff);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        %                         q = varargin{1};s = varargin{2};
                        %                         [dP, dPdQ, dPdS] = DuctNetwork.Fan(q,s,DuctNetwork.Table_NVT250.LnP_Coeff,DuctNetwork.Table_NVT250.LnQ_Coeff,DuctNetwork.Table_NVT250.Eff_max);
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        %                         q = varargin{1};s = varargin{2};
                        %                         [dP, dPdQ, dPdS] = DuctNetwork.Fan(q,s,DuctNetwork.Table_NVT250.LnP_Coeff,DuctNetwork.Table_NVT250.LnQ_Coeff,DuctNetwork.Table_NVT250.Eff_max);
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Fan Model VNT250';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.Fan_NVT250};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'Fan Speed(RPM)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-3];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = FanQuadratic_VSD(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = -s(1)*s(3)^2;
                        else
                            varargout{ii} = s(1)*(q^2/s(2)^2-s(3)^2);
                        end
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = 0;
                        else
                            varargout{ii} = s(1)*2*q/s(2)^2;
                        end
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = [-s(3)^2,0,-2*s(1)*s(3)];
                        else
                            varargout{ii} = [q^2/s(2)^2-s(3)^2, -s(1)*2*q^2/s(2)^3,-2*s(1)*s(3)];
                        end
                    case 'Model_Description'
                        varargout{ii}='Fan Using Quadratic Model and Fan Affinity Law';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.FanQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:3};
                    case 'Parameter_Description'
                        varargout{ii}={'Max Pressure(Pa)','Max Flow(m^3/s)','Scale Factor'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-3,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e3,1e6];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = PressureSource(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        varargout{ii} = -s;
                    case 'dPdQ'
                        varargout{ii} = 0;
                    case 'dPdS'
                        varargout{ii} = -1;
                    case 'Model_Description'
                        varargout{ii}='Flow source with constant pressure';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.PressureSource};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'Pressure(Pa)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-3];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = Louver(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        varargout{ii} = s(1)*tanh(q/3/s(2));
                    case 'dPdQ'
                        varargout{ii} = s(1)*(1-tanh(q/3/s(2))^2)/3/s(2);
                    case 'dPdS'
                        varargout{ii} = [tanh(q/3/s(2)), -s(1)*(1-tanh(q/3/s(2))^2)*q/3/s(2)^2];
                    case 'Model_Description'
                        varargout{ii}='Louver with constant pressure drop';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.Louver};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Pressure(Pa)','Saturated Flow (m^3/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-3, 1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9, 1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false, false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = DuctQuadratic(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = s*q*abs(q);
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = s*2*abs(q);
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = q*abs(q);
                    case 'Model_Description'
                        varargout{ii}='Duct Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.DuctQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'Resistance(Pa/(m^3/s)^2)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e18];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = DuctQuadraticExp(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = exp(s)*q*abs(q);
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = exp(s)*2*abs(q);
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = exp(s)*q*abs(q);
                    case 'Model_Description'
                        varargout{ii}='Duct Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.DuctQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'log of Resistance(Pa/(m^3/s)^2)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[-2000];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[2000];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularFlexible(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:nargout
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        K = s(2);
                        D = s(3);
                        rho = s(4);
                        e = s(5);
                        nu = s(6);
                        [Pdrop_, dPdQ_, dPdS_] = DuctNetwork.CircularDarcyWeisbach({'Pdrop','dPdQ','dPdS'},q,s([1,3:6]));
                        PDCF = 1+0.58*K*exp(-4.96*D);
                        dP = Pdrop_*PDCF;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        dPdq = dPdQ_*PDCF;
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        dPdL = dPdS_(1)*PDCF;
                        dPdK = Pdrop_*PDCF/K;
                        dPdD = dPdS_(2)*PDCF-Pdrop_*(PDCF-1)*4.96;
                        dPdrho = dPdS_(3)*PDCF;
                        dPde = dPdS_(4)*PDCF;
                        dPdnu = dPdS_(5)*PDCF;
                        varargout{ii}=[dPdL, dPdK, dPdD, dPdrho, dPde, dPdnu];
                    case 'Model_Description'
                        varargout{ii}='Circular Flexible Duct';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularFlexible};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Compression Percent(%)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-6,1e-9,1e-9,1e-6,1e-6,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e12,1,1e9,1e4,1e6,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularDarcyWeisbach(query, varargin)
            if ischar(query)
                query = {query};
            end
            varargout = cell(1,nargout);
            [varargout{:}] = DuctNetwork.CircularDarcyWeisbachHaaland(query, varargin{:});
        end
        
        function varargout = CircularDarcyWeisbachHaaland(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:nargout
                switch query{ii}
                    case 'Pdrop'
                        %                         if (~exist('dP','var'))
                        %                             DataInitiation();
                        %                         end;
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        D = s(2);
                        rho = s(3);
                        e = s(4);
                        nu = s(5);
                        
                        %                       Area = pi*(D/2)^2;
                        %                       V = q/Area;
                        %                       Re = abs(V)*D/nu;
                        Re = 4*abs(q)/nu/pi/D;
                        
                        lambda = 1/(1+exp(-(Re-3750)/250));
                        Cf_lam = 64/Re;
                        A = (e/3.7/D)^3.33;
                        B = (6.9/Re)^3;
                        T = log10(A+B);
                        Cf_turb = (-0.6*T)^(-2);
                        Cf = (1-lambda)*Cf_lam + lambda*Cf_turb;
                        % M = L/D/2*rho*V*abs(V);
                        M = 8*L*rho/D^5/pi^2*q*abs(q);
                        dP = Cf*M;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        %                         if (~exist('dP','var'))
                        %                             DataInitiation();
                        %                         end;
                        %                         if (~exist('dCf_turbdAB','var'))
                        %                             dCf_turbdAB = -1/0.18/T^3/log(10)/(A+B);
                        %                         end;
                        dCf_lamdabsq = -Cf_lam/abs(q);
                        dCf_turbdAB = -2*Cf_turb/T/log(10)/(A+B);
                        dCf_turbdabsq = -dCf_turbdAB*3*B/abs(q);
                        dlambdadabsq = lambda*(1-lambda)/250*Re/abs(q);
                        dCfdabsq = (Cf_turb-Cf_lam)*dlambdadabsq + (1-lambda)*dCf_lamdabsq + lambda*dCf_turbdabsq;
                        dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        %                         if (~exist('dP','var'))
                        %                             DataInitiation();
                        %                         end;
                        %                         if (~exist('dCf_turbdAB','var'))
                        %                             dCf_turbdAB = -1/0.18/T^3/log(10)/(A+B);
                        %                         end;
                        dPdL = dP/L;
                        dCf_lamdD = Cf_lam/D;
                        dCf_turbdD = dCf_turbdAB*(-3.33*A+3*B)/D;
                        dlambdadD = -lambda*(1-lambda)/250*Re/D;
                        dCfdD = (Cf_turb-Cf_lam)*dlambdadD + (1-lambda)*dCf_lamdD + lambda*dCf_turbdD;
                        dPdD = M*(dCfdD-5*Cf/D);
                        dPdrho = dP/rho;
                        dCf_turbde = dCf_turbdAB*3.33*A/e;
                        dPde = M*lambda*dCf_turbde;
                        dCf_lamdnu = Cf_lam/nu;
                        dCf_turbdnu = dCf_turbdAB*3*B/nu;
                        dlambdadnu = -lambda*(1-lambda)/250*Re/nu;
                        dCfdnu = (Cf_turb-Cf_lam)*dlambdadnu + (1-lambda)*dCf_lamdnu + lambda*dCf_turbdnu;
                        dPdnu = M*dCfdnu;
                        varargout{ii}=[dPdL, dPdD, dPdrho, dPde, dPdnu];
                    case 'Model_Description'
                        varargout{ii}='Circular Straight Duct Using Darcy Weisbach Equation by Haaland Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-6,1e-9,1e-6,1e-6,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e12,1e9,1e4,1e6,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularDarcyWeisbachChurchill(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        D = s(2);
                        rho = s(3);
                        e = s(4);
                        nu = s(5);
                        
                        %                       Area = pi*(D/2)^2;
                        %                       V = q/Area;
                        %                       Re = abs(V)*D/nu;
                        Re = 4*abs(q)/nu/pi/D;
                        
                        T1 = power((7/Re),0.9);
                        T2 = T1 +(0.27*e/D);
                        T3 = -2.457*log(T2);
                        A = T3^16;
                        B = power((37530/Re),16);
                        T4 = power((8/Re),12);
                        T5 = power(A+B,-1.5);
                        Cf = 8*power(T4+T5,1/12);
                        % M = L/D/2*rho*V*abs(V);
                        M = 8*L*rho/D^5/pi^2*q*abs(q);
                        dP = Cf*M;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        dAdabsq = 35.3808*A/T3/T2*T1/abs(q);
                        dBdabsq = -16*B/abs(q);
                        dT4dabsq = -12*T4/abs(q);
                        dT5dabsq = -1.5*T5/(A+B)*(dAdabsq+dBdabsq);
                        dCfdabsq = 1/12*Cf/(T4+T5)*(dT4dabsq+dT5dabsq);
                        dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        dPdL = dP/L;
                        dT4dD = 12*T4/D;
                        dAdD = -2.457*16*A/T3/T2*(0.9*T1/D-0.27*e/D^2);
                        dBdD = 16*B/D;
                        dT5dD = -1.5*T5/(A+B)*(dAdD+dBdD);
                        dCfdD = 1/12*Cf/(T4+T5)*(dT4dD+dT5dD);
                        dPdD = M*(dCfdD-5*Cf/D);
                        dPdrho = dP/rho;
                        % dAde = -2.457*16*0.27*A/T3/T2/D;
                        % dT5de = -1.5*T5/(A+B)*dAde;
                        % dCfde = 1/12*Cf/(T4+T5)*dT5de;
                        % dPde = M*dCfde;
                        dPde = 1.32678*dP*T3^15/T2/(T4+T5)/(A+B)^2.5/D;
                        dT4dnu = 12*T4/nu;
                        dAdnu = -2.457*16*0.9*A/T3/T2*T1/nu;
                        dBdnu = 16*B/nu;
                        dT5dnu = -1.5*T5/(A+B)*(dAdnu+dBdnu);
                        dCfdnu = 1/12*Cf/(T4+T5)*(dT4dnu+dT5dnu);
                        dPdnu = M*dCfdnu;
                        dPds = [dPdL, dPdD, dPdrho, dPde, dPdnu];
                        varargout{ii}=dPds;
                    case 'Model_Description'
                        varargout{ii}='Circular Straight Duct Using Darcy Weisbach Equation by Churchill Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-6,1e-9,1e-6,1e-6,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e12,1e9,1e4,1e6,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = RectangularDarcyWeisbach(query, varargin)
            if ischar(query)
                query = {query};
            end
            varargout = cell(1,nargout);
            [varargout{:}] = DuctNetwork.RectangularDarcyWeisbachHaaland(query, varargin{:});
        end
        
        function varargout = RectangularDarcyWeisbachHaaland(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:nargout
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        H = s(2);
                        W = s(3);
                        Dh = 1.3*H^0.625*W^0.625/(H+W)^0.25;
                        rho = s(4);
                        e = s(5);
                        nu = s(6);
                        % Re = abs(V)*Dh/nu;
                        Re = abs(q)*Dh/H/W/nu;
                        
                        lambda = 1/(1+exp(-(Re-3750)/250));
                        Cf_lam = 64/Re;
                        A = (e/3.7/Dh)^3.33;
                        B = (6.9/Re)^3;
                        T = log10(A+B);
                        Cf_turb = (-0.6*T)^(-2);
                        Cf = (1-lambda)*Cf_lam + lambda*Cf_turb;
                        % M = L/Dh/2*rho*V*abs(V);
                        M = L*rho*q*abs(q)/Dh/2/H^2/W^2;
                        dP = Cf*M;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        dCf_lamdabsq = -Cf_lam/abs(q);
                        dCf_turbdAB = -2*Cf_turb/T/log(10)/(A+B);
                        dCf_turbdabsq = -dCf_turbdAB*3*B/abs(q);
                        dlambdadabsq = lambda*(1-lambda)/250*Re/abs(q);
                        dCfdabsq = (Cf_turb-Cf_lam)*dlambdadabsq + (1-lambda)*dCf_lamdabsq + lambda*dCf_turbdabsq;
                        dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        dPdL = dP/L;
                        Kh = 0.625/H-0.25/(H+W);
                        Kh2 = Kh-1/H;
                        % dDhdH = Dh*Kh;
                        dRedH  = Re*Kh2;
                        dCf_lamdH = -Cf_lam*Kh2;
                        dAdH = -3.33*A*Kh;
                        dBdH = -3*B*Kh2;
                        dCf_turbdH = dCf_turbdAB*(dAdH+dBdH);
                        dlambdadH = lambda*(1-lambda)/250*dRedH;
                        dCfdH = (Cf_turb-Cf_lam)*dlambdadH + (1-lambda)*dCf_lamdH + lambda*dCf_turbdH;
                        dMdH = M*(2*Kh2-3*Kh);
                        dPdH = dCfdH*M+Cf*dMdH;
                        Kw = 0.625/W-0.25/(H+W);
                        Kw2 = Kw-1/W;
                        % dDhdW = Dh*Kw;
                        dRedW  = Re*Kw2;
                        dCf_lamdW = -Cf_lam*Kw2;
                        dAdW = -3.33*A*Kw;
                        dBdW = -3*B*Kw2;
                        dCf_turbdW = dCf_turbdAB*(dAdW+dBdW);
                        dlambdadW = lambda*(1-lambda)/250*dRedW;
                        dCfdW = (Cf_turb-Cf_lam)*dlambdadW + (1-lambda)*dCf_lamdW + lambda*dCf_turbdW;
                        dMdW = M*(2*Kw2-3*Kw);
                        dPdW = dCfdW*M+Cf*dMdW;
                        dPdrho = dP/rho;
                        dCf_turbde = dCf_turbdAB*3.33*A/e;
                        dPde = M*lambda*dCf_turbde;
                        dCf_lamdnu = Cf_lam/nu;
                        dCf_turbdnu = dCf_turbdAB*3*B/nu;
                        dlambdadnu = -lambda*(1-lambda)/250*Re/nu;
                        dCfdnu = (Cf_turb-Cf_lam)*dlambdadnu + (1-lambda)*dCf_lamdnu + lambda*dCf_turbdnu;
                        dPdnu = M*dCfdnu;
                        varargout{ii}=[dPdL, dPdH, dPdW, dPdrho, dPde, dPdnu];
                    case 'Model_Description'
                        varargout{ii}='Rectangular Straight Duct Using Darcy Weisbach Equation by Haaland Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Height(m)','Width(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-6,1e-9,1e-9,1e-6,1e-6,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e12,1e9,1e9,1e4,1e6,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = RectangularDarcyWeisbachChurchill(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        H = s(2);
                        W = s(3);
                        Dh = 1.3*H^0.625*W^0.625/(H+W)^0.25;
                        rho = s(4);
                        e = s(5);
                        nu = s(6);
                        % Re = abs(V)*Dh/nu;
                        Re = abs(q)*Dh/H/W/nu;
                        
                        T1 = power((7/Re),0.9);
                        T2 = T1 +(0.27*e/Dh);
                        T3 = -2.457*log(T2);
                        A = T3^16;
                        B = power((37530/Re),16);
                        T4 = power((8/Re),12);
                        T5 = power(A+B,-1.5);
                        Cf = 8*power(T4+T5,1/12);
                        % M = L/Dh/2*rho*V*abs(V);
                        M = L*rho*q*abs(q)/Dh/2/H^2/W^2;
                        dP = Cf*M;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        dAdabsq = 35.3808*A/T3/T2*T1/abs(q);
                        dBdabsq = -16*B/abs(q);
                        dT4dabsq = -12*T4/abs(q);
                        dT5dabsq = -1.5*T5/(A+B)*(dAdabsq+dBdabsq);
                        dCfdabsq = 1/12*Cf/(T4+T5)*(dT4dabsq+dT5dabsq);
                        dPdq = M*(dCfdabsq*sign(q)+2*Cf/q);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        dPdL = dP/L;
                        Kh = 0.625/H-0.25/(H+W);
                        Kh2 = Kh-1/H;
                        % dDhdH = Dh*Kh;
                        % dRedH  = Re*Kh2;
                        dT4dH = -12*T4*Kh2;
                        dT2dH = -0.9*T1*Kh2+0.27*e/Dh*(0.25/(H+W)-0.625/H);
                        dAdH = -2.457*16*A/T3/T2*dT2dH;
                        dBdH = -16*B*Kh2;
                        dT5dH = -1.5*T5/(A+B)*(dAdH+dBdH);
                        dCfdH = 1/12*Cf/(T4+T5)*(dT4dH+dT5dH);
                        dMdH = M*(0.25/(H+W)-2.625/H);
                        dPdH = dCfdH*M+Cf*dMdH;
                        Kw = 0.625/W-0.25/(H+W);
                        Kw2 = Kw-1/W;
                        % dDhdW = Dh*Kw;
                        % dRedW  = Re*Kw2;
                        dT4dW = -12*T4*Kw2;
                        dT2dW = -0.9*T1*Kw2+0.27*e/Dh*(0.25/(H+W)-0.625/W);
                        dAdW = -2.457*16*A/T3/T2*dT2dW;
                        dBdW = -16*B*Kw2;
                        dT5dW = -1.5*T5/(A+B)*(dAdW+dBdW);
                        dCfdW = 1/12*Cf/(T4+T5)*(dT4dW+dT5dW);
                        dMdW = M*(0.25/(H+W)-2.625/W);
                        dPdW = dCfdW*M+Cf*dMdW;
                        dPdrho = dP/rho;
                        % dAde = -2.457*16*0.27*A/T3/T2/Dh;
                        % dT5de = -1.5*T5/(A+B)*dAde;
                        % dCfde = 1/12*Cf/(T4+T5)*dT5de;
                        % dPde = M*dCfde;
                        dPde = 1.32678*dP*T3^15/T2/(T4+T5)/(A+B)^2.5/Dh;
                        dT4dnu = 12*T4/nu;
                        dAdnu = -2.457*16*0.9*A/T3/T2*T1/nu;
                        dBdnu = 16*B/nu;
                        dT5dnu = -1.5*T5/(A+B)*(dAdnu+dBdnu);
                        dCfdnu = 1/12*Cf/(T4+T5)*(dT4dnu+dT5dnu);
                        dPdnu = M*dCfdnu;
                        varargout{ii}=[dPdL, dPdH, dPdW, dPdrho, dPde, dPdnu];
                    case 'Model_Description'
                        varargout{ii}='Rectangular Straight Duct Using Darcy Weisbach Equation by Churchill Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Height(m)','Width(m)','Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-6,1e-9,1e-9,1e-6,1e-6,1e-9];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e12,1e9,1e9,1e4,1e6,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularTJunction(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(4,1);
                    case 'Model_Description'
                        varargout{ii}='Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Horizontal,@DuctNetwork.CircularTJunction_Horizontal,@DuctNetwork.CircularTJunction_Vertical};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4];[2,1,3,4];[3,1,2,4]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
        end
        
        function varargout = CircularTJunction_Horizontal(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Horizontal part of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SD5-18, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.SD5_18(abs(q([1,2,3])), s([1,2,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.SD5_9(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ED5-3 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SD5-9 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ED5-3 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.ED5_3(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ED5-4, use Cb
                        if s(1)>=s(2) % D1>D2, use Cb1
                            [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.ED5_4(abs(q([1,2,3])), s([1,2,3,4]), 'b1');
                        else % D1<D2, use Cb2
                            [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_4(abs(q([2,1,3])), s([2,1,3,4]), 'b2');
                        end
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularTJunction_Vertical(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Vertical part of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Vertical};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Diameter(m)','Main 1 Diameter(m)','Main 2 Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.SD5_9(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.SD5_9(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ED5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SD5-18 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.ED5_3(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the main 1, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_3(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularY30Junction(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(4,1);
                    case 'Model_Description'
                        varargout{ii}='Circular Y-30-Junction using ED5-1,ED5-4,SD5-3,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY30Junction_Horizontal,@DuctNetwork.CircularY30Junction_Horizontal,@DuctNetwork.CircularY30Junction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4];[2,1,3,4];[3,1,2,4]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
        end
        
        function varargout = CircularY30Junction_Horizontal(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Horizontal part of Circular Y-30-Junction using ED5-1,ED5-4,SD5-3,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY30Junction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SD5-18, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.SD5_18(abs(q([1,2,3])), s([1,2,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-3 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.SD5_3(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ED5-1 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SD5-3 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ED5-1 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.ED5_1(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ED5-4, use Cb
                        if s(1)>=s(2) % D1>D2, use Cb1
                            [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.ED5_4(abs(q([1,2,3])), s([1,2,3,4]), 'b1');
                        else % D1<D2, use Cb2
                            [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_4(abs(q([2,1,3])), s([2,1,3,4]), 'b2');
                        end
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularY30Junction_Side(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Side branch of Circular Y-30-Junction using ED5-1,ED5-4,SD5-3,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY30Junction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Diameter(m)','Main 1 Diameter(m)','Main 2 Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SD5-3 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.SD5_3(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-3 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.SD5_3(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ED5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SD5-18 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ED5-1 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.ED5_1(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, T converge ED5-1 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_1(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularY45Junction(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(4,1);
                    case 'Model_Description'
                        varargout{ii}='Circular Y-45-Junction using ED5-2,ED5-4,SD5-1,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY45Junction_Horizontal,@DuctNetwork.CircularY45Junction_Horizontal,@DuctNetwork.CircularY45Junction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4];[2,1,3,4];[3,1,2,4]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
        end
        
        function varargout = CircularY45Junction_Horizontal(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Horizontal part of Circular Y-45-Junction using ED5-2,ED5-4,SD5-1,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY45Junction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SD5-18, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.SD5_18(abs(q([1,2,3])), s([1,2,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-1 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.SD5_1(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ED5-2 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SD5-1 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ED5-2 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.ED5_2(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ED5-4, use Cb
                        if s(1)>=s(2) % D1>D2, use Cb1
                            [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.ED5_4(abs(q([1,2,3])), s([1,2,3,4]), 'b1');
                        else % D1<D2, use Cb2
                            [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_4(abs(q([2,1,3])), s([2,1,3,4]), 'b2');
                        end
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularY45Junction_Side(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Side branch of Circular Y-45-Junction using ED5-2,ED5-4,SD5-1,SD5-18';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularY45Junction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Diameter(m)','Main 1 Diameter(m)','Main 2 Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SD5-1 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.SD5_1(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-1 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.SD5_1(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ED5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SD5-18 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ED5-2 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.ED5_2(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, T converge ED5-2 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.ED5_2(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = RectangularTJunction(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(7,1);
                    case 'Model_Description'
                        varargout{ii}='Rectangular T-Junction using SR5-13,SR5-15,ER5-3,ER5-5';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularTJunction_Horizontal,@DuctNetwork.RectangularTJunction_Horizontal,@DuctNetwork.RectangularTJunction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Height (m)','Main 1 Width(m)','Main 2 Height (m)','Main 2 Width(m)','Branch Height (m)','Branch Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e9,1e9,1e4];
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4,5,6,7];[3,4,1,2,5,6,7];[5,6,1,2,3,4,7]};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0,0];
                end
            end
        end
        
        function varargout = RectangularTJunction_Horizontal(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Horizontal part of Rectangular T-Junction using SR5-13,SR5-15,ER5-3,ER5-5';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularTJunction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:7};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Height (m)','Main 1 Width(m)','Main 2 Height (m)','Main 2 Width(m)','Branch Height (m)','Branch Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SR5-15, use Cb
                        % H = (H1+H2+H3)/3; Wb1 = W1; Wb2 = W2; Wc = W3; [H;Wb1;Wb2;Wc;rho]
                        Convertion = sparse([1,1,1,2,3,4,5],[1,3,5,2,4,6,7],[1/3,1/3,1/3,1,1,1,1],5,7);
                        [dP, dPdQ([1,2,3]), dPdS] = DuctNetwork.SR5_15(abs(q([1,2,3])),s*Convertion','b');
                        dPdS = dPdS*Convertion;
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SR5-13 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,2,5,6,3,4,7])] = DuctNetwork.SR5_13(abs(q([1,3,2])),s([1,2,5,6,3,4,7]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ER5-3 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SR5-13 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ER5-3 at upstream side, use Cs
                        % H = (H1+H2+H3)/3; Ws=W1; Wb=W3; Wc = W2; [H;Ws;Wb;Wc;rho]
                        Convertion = sparse([1,1,2,2,3,4,5],[1,3,2,4,5,6,7],[1/2,1/2,1/2,1/2,1,1,1],5,7);
                        [dP, dPdQ([1,3,2]), dPdS] = DuctNetwork.ER5_3(abs(q([1,3,2])),s*Convertion','s');
                        dPdS = dPdS*Convertion;
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ER5-5, use Cb
                        % H = (H1+H2+H3)/3; Wb1 = W1; Wb2 = W2; Wc = W3; [H;Wb1;Wb2;Wc;rho]
                        Convertion = sparse([1,1,1,2,3,4,5],[1,3,5,2,4,6,7],[1/3,1/3,1/3,1,1,1,1],5,7);
                        [dP, dPdQ([1,2,3]), dPdS] = DuctNetwork.ER5_5(abs(q([1,2,3])),s*Convertion','b');
                        dPdS = dPdS*Convertion;
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = RectangularTJunction_Side(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Side branch of Rectangular T-Junction using SR5-13,SR5-15,ER5-3,ER5-5';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularTJunction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:7};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Height (m)','Branch Width(m)','Main 1 Height (m)','Main 1 Width(m)','Main 2 Height (m)','Main 2 Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SR5-13 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([3,4,1,2,5,6,7])] = DuctNetwork.SR5_13(abs(q([2,1,3])),s([3,4,1,2,5,6,7]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SR5-13 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([5,6,1,2,3,4,7])] = DuctNetwork.SR5_13(abs(q([3,1,2])),s([5,6,1,2,3,4,7]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ER5-5 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SR5-15 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ER5-3 at branch side, use Cb
                        % H = (H1+H2+H3)/3; Ws = W3; Wb = W1; Wc = W2; [H;Ws;Wb;Wc;rho]
                        Convertion = sparse([1,1,2,2,3,4,5],[1,3,2,4,5,6,7],[1/2,1/2,1/2,1/2,1,1,1],5,7);
                        [dP, dPdQ([3,1,2]), dPdS] = DuctNetwork.ER5_3(abs(q([3,1,2])),s*Convertion','b');
                        dPdS = dPdS*Convertion;
                    case 6 %[1,1,0]*[4;2;1], + + -, T converge ER5-3 at branch side, use Cb
                        % H = (H1+H2+H3)/3; Ws = W2; Wb = W1; Wc = W3; [H;Ws;Wb;Wc;rho]
                        Convertion = sparse([1,1,2,2,3,4,5],[1,3,2,4,5,6,7],[1/2,1/2,1/2,1/2,1,1,1],5,7);
                        [dP, dPdQ([2,1,3]), dPdS] = DuctNetwork.ER5_3(abs(q([2,1,3])),s*Convertion','b');
                        dPdS = dPdS*Convertion;
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = RectangularYJunction(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(5,1);
                    case 'Model_Description'
                        varargout{ii}='Rectangular Y-Junction using SR5-1,ER5-1,SR5-14,ER5-4';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularYJunction_Horizontal,@DuctNetwork.RectangularYJunction_Horizontal,@DuctNetwork.RectangularYJunction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4,5];[1,3,2,4,5];[1,4,2,3,5]};
                    case 'Parameter_Description'
                        varargout{ii}={'Height (m)','Main 1 Width(m)','Main 2 Width(m)','Branch Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                end
            end
        end
        
        function varargout = RectangularYJunction_Horizontal(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Horizontal part of Rectangular Y-Junction using SR5-1,ER5-1,SR5-14,ER5-4';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularYJunction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height (m)','Main 1 Width(m)','Main 2 Width(m)','Branch Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SR5-14, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4,5])] = DuctNetwork.SR5_14(abs(q([1,2,3])), s([1,2,3,4,5]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SR5-1 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,2,4,3,5])] = DuctNetwork.SR5_1(abs(q([1,3,2])),s([1,2,4,3,5]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ER5-1 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SR5-1 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ER5-1 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,2,4,3,5])] = DuctNetwork.ER5_1(abs(q([1,3,2])),s([1,2,4,3,5]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ER5-4, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4,5])] = DuctNetwork.ER5_4(abs(q([1,2,3])), s([1,2,3,4,5]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = RectangularYJunction_Side(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = Calculation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Side branch of Rectangular Y-Junction using SR5-1,ER5-1,SR5-14,ER5-4';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.RectangularYJunction_Side};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height (m)','Main 1 Width(m)','Main 2 Width(m)','Branch Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=Calculation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SR5-1 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([1,3,2,4,5])] = DuctNetwork.SR5_1(abs(q([2,1,3])),s([1,3,2,4,5]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SR5-1 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([1,4,2,3,5])] = DuctNetwork.SR5_1(abs(q([3,1,2])),s([1,4,2,3,5]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ER5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SR5-14 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ER5-1 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([1,4,2,3,5])] = DuctNetwork.ER5_1(abs(q([3,1,2])),s([1,4,2,3,5]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, T converge ER5-1 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([1,3,2,4,5])] = DuctNetwork.ER5_1(abs(q([2,1,3])),s([1,3,2,4,5]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CD3_6(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             rho = s(2);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_CD3_6.D};
                        Co_Table = DuctNetwork.Table_CD3_6.Co;
                        gExp = 0.5*s(2)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(2)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(2)*q*abs(q)/s(1)^5/pi^2,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = s(1);
                        dZdq = 0;
                        dZds = [1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CD3_6 Elbow, Pleated, 60 Degree';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD3_6};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CD3_9(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             rho = s(2);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_CD3_9.D};
                        Co_Table = DuctNetwork.Table_CD3_9.Co;
                        gExp = 0.5*s(2)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(2)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(2)*q*abs(q)/s(1)^5/pi^2,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = s(1);
                        dZdq = 0;
                        dZds = [1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CD3_9 Elbow, 5 Gore, 90 Degree';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD3_9};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CD3_17(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             rho = s(2);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_CD3_17.D};
                        Co_Table = DuctNetwork.Table_CD3_17.Co;
                        gExp = 0.5*s(2)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(2)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(2)*q*abs(q)/s(1)^5/pi^2,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = [s(1)];
                        dZdq = [0];
                        dZds = [1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CD3_17 Elbow, Mitered, 45 Degree';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD3_17};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CD6_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             Do = s(1);
                        %             A1/Ao = s(2);
                        %             n = s(3);
                        %             rho = s(4);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_CD6_1.n,DuctNetwork.Table_CD6_1.A1Ao};
                        Co_Table = DuctNetwork.Table_CD6_1.Co;
                        gExp = 0.5*s(4)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(4)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(4)*q*abs(q)/s(1)^5/pi^2,0,0,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = [s(3);s(2)];
                        dZdq = [0;0];
                        dZds = [0,0,1,0;0,1,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CD6_1 Screen';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD6_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Area Ratio(1)','Free Area Ratio(1)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,0,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e3,1e3,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CD9_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             Do = s(1);
                        %             D = s(2);
                        %             Theta = s(3);
                        %             rho = s(4);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_CD9_1.Theta,DuctNetwork.Table_CD9_1.DDo};
                        Co_Table = DuctNetwork.Table_CD9_1.Co;
                        gExp = 0.5*s(4)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(4)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(4)*q*abs(q)/s(1)^5/pi^2,0,0,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = [s(3);s(2)/s(1)];
                        dZdq = [0;0];
                        dZds = [0,0,1,0;-s(2)/s(1)^2,1/s(1),0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %                         if (s(3)<5)*(s(3)>0)+(s(3)<0)+(s(3)>89)*(s(3)<90)+(s(3)>90)
                        %                             disp([s(3),dPdS(3)])
                        %                         end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CD9_1 Damper, Butterfly';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD9_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Duct Diameter(m)','Plate Diameter(m)','Angle(deg)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,90,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,1,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CD9_3(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             rho = s(2);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        varargout{ii}=0.12*0.5*s(2)*q*abs(q)/pi^2/(s(1)/2)^4;
                    case 'dPdQ'
                        varargout{ii}=0.12*s(2)*abs(q)/pi^2/(s(1)/2)^4;
                    case 'dPdS'
                        varargout{ii}=[-0.12*32*s(2)*q*abs(q)/pi^2/s(1)^5,0.12*0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                    case 'Model_Description'
                        varargout{ii}='CD9_3 Fire Damper';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CD9_3};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR3_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             r/W = s(3);
                        %             Theta = s(4);
                        %             rho = s(5);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        GridVector1 = {DuctNetwork.Table_CR3_1.HW,DuctNetwork.Table_CR3_1.rW};
                        InterpTable1 = DuctNetwork.Table_CR3_1.Cp;
                        Z1Exp = [s(1)/s(2);s(3)];
                        dZ1dq = [0;0];
                        dZ1ds = [1/s(2),-s(1)/s(2)^2,0,0,0;0,0,1,0,0];
                        GridVector2 = {DuctNetwork.Table_CR3_1.Theta};
                        InterpTable2 = DuctNetwork.Table_CR3_1.K;
                        Z2Exp = s(4);
                        dZ2dq = 0;
                        dZ2ds = [0,0,0,1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Double_Interp_Gradient(GridVector1, InterpTable1, GridVector2, InterpTable2, Z1Exp, dZ1dq, dZ1ds, Z2Exp, dZ2dq, dZ2ds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR3_1 Elbow, Smooth Radius';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR3_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Radius Ratio(1)','Angle(deg)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-3,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e3,90,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR3_3(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             r/W = s(3);
                        %             Theta = s(4);
                        %             rho = s(5);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        GridVector1 = {DuctNetwork.Table_CR3_3.HW,DuctNetwork.Table_CR3_3.rW};
                        InterpTable1 = DuctNetwork.Table_CR3_3.Cp;
                        Z1Exp = [s(1)/s(2);s(3)];
                        dZ1dq = [0;0];
                        dZ1ds = [1/s(2),-s(1)/s(2)^2,0,0,0;0,0,1,0,0];
                        GridVector2 = {DuctNetwork.Table_CR3_3.Theta};
                        InterpTable2 = DuctNetwork.Table_CR3_3.K;
                        Z2Exp = s(4);
                        dZ2dq = 0;
                        dZ2ds = [0,0,0,1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Double_Interp_Gradient(GridVector1, InterpTable1, GridVector2, InterpTable2, Z1Exp, dZ1dq, dZ1ds, Z2Exp, dZ2dq, dZ2ds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR3_3 Elbow, Smooth Radius';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR3_3};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Radius Ratio(1)','Angle(deg)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-3,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e3,90,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR3_6(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             Theta = s(3);
                        %             rho = s(4);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_CR3_6.HW,DuctNetwork.Table_CR3_6.Theta};
                        Co_Table = DuctNetwork.Table_CR3_6.Co;
                        gExp = 0.5*s(4)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(4)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(4)*q*abs(q)/s(1)^3/s(2)^2,-s(4)*q*abs(q)/s(1)^2/s(2)^3,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        ZExp = [s(1)/s(2);s(3)];
                        dZdq = [0;0];
                        dZds = [1/s(2),-s(1)/s(2)^2,0,0;0,0,1,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR3_6 Damper, Elbow, Mitered';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR3_6};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Angle(deg)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,90,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR3_10(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             rho = s(3);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        varargout{ii}=0.12*0.5*s(3)*q*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdQ'
                        varargout{ii}=0.12*s(3)*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdS'
                        varargout{ii}=[-0.12*s(3)*q*abs(q)/s(1)^3/s(2)^2,-0.12*s(3)*q*abs(q)/s(1)^2/s(2)^3,0.12*0.5*q*abs(q)/s(1)^2/s(2)^2];
                    case 'Model_Description'
                        varargout{ii}='CR3_10 Elbow, Mitered, 90 Degree';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR3_10};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:3};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR3_17(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             L = s(3);
                        %             rho = s(4);
                        %             nu = s(5);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        gExp = 0.5*s(4)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(4)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(4)*q*abs(q)/s(1)^3/s(2)^2,-s(4)*q*abs(q)/s(1)^2/s(2)^3,0,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                        GridVector1 = {DuctNetwork.Table_CR3_17.LW,DuctNetwork.Table_CR3_17.HW};
                        InterpTable1 = DuctNetwork.Table_CR3_17.Cp;
                        Z1Exp = [s(3)/s(2);s(1)/s(2)];
                        dZ1dq = [0;0];
                        dZ1ds = [0,-s(3)/s(2)^2,1/s(2),0,0;1/s(2),-s(1)/s(2)^2,0,0,0];
                        GridVector2 = {DuctNetwork.Table_CR3_17.Re};
                        InterpTable2 = DuctNetwork.Table_CR3_17.Kr;
                        Z2Exp = 2*abs(q)/(s(1)+s(2))/s(5);
                        dZ2dq = 2*sign(q)/(s(1)+s(2))/s(5);
                        dZ2ds = [-2*abs(q)/(s(1)+s(2))^2/s(5),-2*abs(q)/(s(1)+s(2))^2/s(5),0,0,-2*abs(q)/(s(1)+s(2))/s(5)^2];
                        [dP, dPdQ, dPdS] = DuctNetwork.Double_Interp_Gradient(GridVector1, InterpTable1, GridVector2, InterpTable2, Z1Exp, dZ1dq, dZ1ds, Z2Exp, dZ2dq, dZ2ds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR3_17 Elbow, Z Shape';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR3_17};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Length(m)','Density(kg/m^3)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6,1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR6_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             A1/Ao = s(3);
                        %             n = s(4);
                        %             rho = s(5);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_CR6_1.n,DuctNetwork.Table_CR6_1.A1Ao};
                        Co_Table = DuctNetwork.Table_CR6_1.Co;
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        ZExp = [s(4);s(3)];
                        dZdq = [0;0];
                        dZds = [0,0,0,1,0;0,0,1,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR6_1 Screen';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR6_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Area Ratio(1)','Free Area Ratio(1)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e3,1e3,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR6_4(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             d = s(3);
                        %             y = s(4);
                        %             rho = s(5);
                        %             nu = s(6);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_CR6_4.SmAo,DuctNetwork.Table_CR6_4.Re,DuctNetwork.Table_CR6_4.yH};
                        Co_Table = DuctNetwork.Table_CR6_4.Co;
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                        ZExp = [s(3)/s(1);abs(q)*s(3)/s(1)/s(2)/s(6);s(4)/s(1)];
                        dZdq = [0;sign(q)*s(3)/s(1)/s(2)/s(6);0];
                        dZds = [-s(3)/s(1)^2,0,1/s(1),0,0,0;-abs(q)*s(3)/s(1)^2/s(2)/s(6),-abs(q)*s(3)/s(1)/s(2)^2/s(6),abs(q)/s(1)/s(2)/s(6),0,0,-abs(q)*s(3)/s(1)/s(2)/s(6)^2;-s(4)/s(1)^2,0,0,1/s(1),0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR6_4 Obstruction';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR6_4};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Diameter(m)','Distance(m)','Density(kg/m^3)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6,1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR9_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             Theta = s(3);
                        %             rho = s(4);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_CR9_1.Theta,DuctNetwork.Table_CR9_1.HW};
                        Co_Table = DuctNetwork.Table_CR9_1.Co;
                        gExp = 0.5*s(4)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(4)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(4)*q*abs(q)/s(1)^3/s(2)^2,-s(4)*q*abs(q)/s(1)^2/s(2)^3,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        ZExp = [s(3);s(1)/s(2)];
                        dZdq = [0;0];
                        dZds = [0,0,1,0;1/s(2),-s(1)/s(2)^2,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %                         if (s(3)<5)*(s(3)>0)+(s(3)<0)+(s(3)>89)*(s(3)<90)+(s(3)>90)
                        %                             disp([s(3),dPdS(3)])
                        %                         end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='CR9_1 Damper, Butterfly';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR9_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Angle(deg)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,90,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,1,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR9_4(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             rho = s(3);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        varargout{ii}=0.18*0.5*s(3)*q*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdQ'
                        varargout{ii}=0.18*s(3)*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdS'
                        varargout{ii}=[-0.18*s(3)*q*abs(q)/s(1)^3/s(2)^2,-0.18*s(3)*q*abs(q)/s(1)^2/s(2)^3,0.18*0.5*q*abs(q)/s(1)^2/s(2)^2];
                    case 'Model_Description'
                        varargout{ii}='CR9_4 Damper, Parallel & Opposed Airfoil Blades, Open';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR9_4};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:3};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CR9_6(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             rho = s(3);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        varargout{ii}=0.19*0.5*s(3)*q*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdQ'
                        varargout{ii}=0.19*s(3)*abs(q)/s(1)^2/s(2)^2;
                    case 'dPdS'
                        varargout{ii}=[-0.19*s(3)*q*abs(q)/s(1)^3/s(2)^2,-0.19*s(3)*q*abs(q)/s(1)^2/s(2)^3,0.19*0.5*q*abs(q)/s(1)^2/s(2)^2];
                    case 'Model_Description'
                        varargout{ii}='CR9_6 Elbow, Fire Damper, Curtain Type';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CR9_6};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:3};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = ED1_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             t = s(2);
                        %             L = s(3);
                        %             rho = s(4);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_ED1_1.LD,DuctNetwork.Table_ED1_1.tD};
                        Co_Table = DuctNetwork.Table_ED1_1.Co;
                        gExp = 0.5*s(4)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(4)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(4)*q*abs(q)/s(1)^5/pi^2,0,0,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = [s(3)/s(1);s(2)/s(1)];
                        dZdq = [0;0];
                        dZds = [-s(3)/s(1)^2,0,1/s(1),0;-s(2)/s(1)^2,1/s(1),0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='ED1_1 Duct mounted in wall';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.ED1_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Wall Thickness(m)','Extension(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e3,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = ED1_3(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             r = s(2);
                        %             rho = s(3);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_ED1_3.rD};
                        Co_Table = DuctNetwork.Table_ED1_3.Co;
                        gExp = 0.5*s(3)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(3)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(3)*q*abs(q)/s(1)^5/pi^2,0,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = s(2)/s(1);
                        dZdq = 0;
                        dZds = [-s(2)/s(1)^2,1/s(1),0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='ED1_3 Bellmouth with wall';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.ED1_3};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:3};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Radius(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = ED7_2(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             D = s(1);
                        %             rD = s(2);
                        %             L = s(3);
                        %             rho = s(4);
                        %             Area = pi*(s(1)/2)^2;
                        %             V = q/pi/(s(1)/2)^2;
                        GridVec = {DuctNetwork.Table_ED7_2.LD,DuctNetwork.Table_ED7_2.rD};
                        Co_Table = DuctNetwork.Table_ED7_2.Co;
                        gExp = 0.5*s(4)*q*abs(q)/pi^2/(s(1)/2)^4;
                        dgdq = s(4)*abs(q)/pi^2/(s(1)/2)^4;
                        dgds = [-32*s(4)*q*abs(q)/s(1)^5/pi^2,0,0,0.5*q*abs(q)/pi^2/(s(1)/2)^4];
                        ZExp = [s(3)/s(1);s(2)];
                        dZdq = [0;0];
                        dZds = [-s(3)/s(1)^2,0,1/s(1),0;0,1,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='ED7_2 Fan Inlet, Centrifugal';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.ED7_2};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Diameter(m)','Radius/Diameter Ratio(1)','Length(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-6,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e6,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = ER4_3(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             D = s(3);
                        %             L = s(4);
                        %             rho = s(5);
                        %             Ao = s(1)*s(2);
                        %             A1 = pi*(s(3)/2)^2;
                        %             Vo = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_ER4_3.tanTheta,DuctNetwork.Table_ER4_3.AoA1};
                        Co_Table = DuctNetwork.Table_ER4_3.Co;
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        if abs(s(3)-s(1)) > abs(s(3)-s(2))% use abs(s(3)-s(1)) in Z expression
                            if s(3) > s(1)
                                ZExp = [(s(3)-s(1))/2/s(4);s(1)*s(2)/pi/(s(3)/2)^2];
                                dZdq = [0;0];
                                dZds = [-1/2/s(4),0,1/2/s(4),-(s(3)-s(1))/2/s(4)^2,0;s(2)/pi/(s(3)/2)^2,s(1)/pi/(s(3)/2)^2,-8*s(1)*s(2)/pi/s(3)^3,0,0];
                            else
                                ZExp = [(s(1)-s(3))/2/s(4);s(1)*s(2)/pi/(s(3)/2)^2];
                                dZdq = [0;0];
                                dZds = [1/2/s(4),0,-1/2/s(4),-(s(1)-s(3))/2/s(4)^2,0;s(2)/pi/(s(3)/2)^2,s(1)/pi/(s(3)/2)^2,-8*s(1)*s(2)/pi/s(3)^3,0,0];
                            end
                        else% use abs(s(3)-s(2)) in Z expression
                            if s(3) > s(2)
                                ZExp = [(s(3)-s(2))/2/s(4);s(1)*s(2)/pi/(s(3)/2)^2];
                                dZdq = [0;0];
                                dZds = [0,-1/2/s(4),1/2/s(4),-(s(3)-s(2))/2/s(4)^2,0;s(2)/pi/(s(3)/2)^2,s(1)/pi/(s(3)/2)^2,-8*s(1)*s(2)/pi/s(3)^3,0,0];
                            else
                                ZExp = [(s(2)-s(3))/2/s(4);s(1)*s(2)/pi/(s(3)/2)^2];
                                dZdq = [0;0];
                                dZds = [0,1/2/s(4),-1/2/s(4),-(s(2)-s(3))/2/s(4)^2,0;s(2)/pi/(s(3)/2)^2,s(1)/pi/(s(3)/2)^2,-8*s(1)*s(2)/pi/s(3)^3,0,0];
                            end
                        end
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='ER4_3 Transition, Rectangular to Round';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.ER4_3};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Diameter(m)','Length(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR2_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             rho = s(3);
                        %             nu = s(4);
                        %             A = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        if 2*abs(q)/(s(1)+s(2))/s(4) < 3000 % laminar flow
                            GridVec = {DuctNetwork.Table_SR2_1.HW};
                            Co_Table = DuctNetwork.Table_SR2_1.Co;
                            gExp = 0.5*s(3)*q*abs(q)/s(1)^2/s(2)^2;
                            dgdq = s(3)*abs(q)/s(1)^2/s(2)^2;
                            dgds = [-s(3)*q*abs(q)/s(1)^3/s(2)^2,-s(3)*q*abs(q)/s(1)^2/s(2)^3,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                            ZExp = [s(1)/s(2)];
                            dZdq = 0;
                            dZds = [1/s(2),-s(1)/s(2)^2,0,0];
                            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        else % turbulent flow, Co=1
                            dP = 0.5*s(3)*q*abs(q)/s(1)^2/s(2)^2;
                            dPdQ = s(3)*abs(q)/s(1)^2/s(2)^2;
                            dPdS = [-s(3)*q*abs(q)/s(1)^3/s(2)^2,-s(3)*q*abs(q)/s(1)^2/s(2)^3,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                        end
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR2_1 Abrupt Exit';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR2_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Density(kg/m^3)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-6,1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e4,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR2_3(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             Ho = s(1);
                        %             Wo = s(2);
                        %             Theta = s(3);
                        %             L = s(4);
                        %             rho = s(5);
                        %             nu = s(6);
                        %             A1/Ao = H1*Wo/(Ho*Wo) = H1/Ho;
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR2_3.Theta,DuctNetwork.Table_SR2_3.Re,DuctNetwork.Table_SR2_3.A1Ao};
                        Co_Table = DuctNetwork.Table_SR2_3.Co;
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                        ZExp = [s(3);2*abs(q)/(s(1)+s(2))/s(6);2*tan(s(3)*pi/360)*s(4)/s(1)+1];
                        dZdq = [0;2*sign(q)/(s(1)+s(2))/s(6);0];
                        dZds = [0,0,1,0,0,0;-2*abs(q)/(s(1)+s(2))^2/s(6),-2*abs(q)/(s(1)+s(2))^2/s(6),0,0,0,-2*abs(q)/(s(1)+s(2))/s(6)^2;-2*tan(s(3)*pi/360)*s(4)/s(1)^2,0,2*s(4)/s(1)*sec(s(3)*pi/360)^2*pi/360,2*tan(s(3)*pi/360)/s(1),0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR2_3 Plain Diffuser, Free Discharge';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR2_3};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Angle(deg)','Length(m)','Density(kg/m^3)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,0,1e-9,1e-6,1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,90,1e9,1e4,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR2_5(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             Ho = s(1);
                        %             Wo = s(2);
                        %             H1 = s(3);
                        %             W1 = s(4);
                        %             L = s(5);
                        %             rho = s(6);
                        %             nu = s(7);
                        %             Ao = s(1)*s(2);
                        %             A1 = s(3)*s(4);
                        %             Vo = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR2_5.tanTheta,DuctNetwork.Table_SR2_5.Re,DuctNetwork.Table_SR2_5.A1Ao};
                        Co_Table = DuctNetwork.Table_SR2_5.Co;
                        gExp = 0.5*s(6)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(6)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(6)*q*abs(q)/s(1)^3/s(2)^2,-s(6)*q*abs(q)/s(1)^2/s(2)^3,0,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2,0];
                        if abs(s(3)-s(1)) > abs(s(4)-s(2))% use abs(s(3)-s(1)) in Z expression
                            if s(3) > s(1)
                                ZExp = [(s(3)-s(1))/2/s(5);2*abs(q)/(s(1)+s(2))/s(7);s(3)*s(4)/s(1)/s(2)];
                                dZdq = [0;2*sign(q)/(s(1)+s(2))/s(7);0];
                                dZds = [-1/2/s(5),0,1/2/s(5),0,-(s(3)-s(1))/2/s(5)^2,0,0;-2*abs(q)/(s(1)+s(2))^2/s(7),-2*abs(q)/(s(1)+s(2))^2/s(7),0,0,0,0,-2*abs(q)/(s(1)+s(2))/s(7)^2;-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0,0];
                            else
                                ZExp = [(s(1)-s(3))/2/s(5);2*abs(q)/(s(1)+s(2))/s(7);s(3)*s(4)/s(1)/s(2)];
                                dZdq = [0;2*sign(q)/(s(1)+s(2))/s(7);0];
                                dZds = [1/2/s(5),0,-1/2/s(5),0,-(s(1)-s(3))/2/s(5)^2,0,0;-2*abs(q)/(s(1)+s(2))^2/s(7),-2*abs(q)/(s(1)+s(2))^2/s(7),0,0,0,0,-2*abs(q)/(s(1)+s(2))/s(7)^2;-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0,0];
                            end
                        else% use abs(s(4)-s(2)) in Z expression
                            if s(4) > s(2)
                                ZExp = [(s(4)-s(2))/2/s(5);2*abs(q)/(s(1)+s(2))/s(7);s(3)*s(4)/s(1)/s(2)];
                                dZdq = [0;2*sign(q)/(s(1)+s(2))/s(7);0];
                                dZds = [0,-1/2/s(5),0,1/2/s(5),-(s(4)-s(2))/2/s(5)^2,0,0;-2*abs(q)/(s(1)+s(2))^2/s(7),-2*abs(q)/(s(1)+s(2))^2/s(7),0,0,0,0,-2*abs(q)/(s(1)+s(2))/s(7)^2;-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0,0];
                            else
                                ZExp = [(s(2)-s(4))/2/s(5);2*abs(q)/(s(1)+s(2))/s(7);s(3)*s(4)/s(1)/s(2)];
                                dZdq = [0;2*sign(q)/(s(1)+s(2))/s(7);0];
                                dZds = [0,1/2/s(5),0,-1/2/s(5),-(s(2)-s(4))/2/s(5)^2,0,0;-2*abs(q)/(s(1)+s(2))^2/s(7),-2*abs(q)/(s(1)+s(2))^2/s(7),0,0,0,0,-2*abs(q)/(s(1)+s(2))/s(7)^2;-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0,0];
                            end
                        end
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR2_5 Pyramidal Diffuser, Free Discharge';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR2_5};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:7};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Height of Discharge(m)','Width of Discharge(m)','Length(m)','Density(kg/m^3)','Kinematic Viscosity(m^2/s)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-9,1e-6,1e-15];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e9,1e4,1e3];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,false,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR2_6(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             W = s(2);
                        %             L = s(3);
                        %             rho = s(4);
                        %             Area = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR2_6.LDh};
                        Co_Table = DuctNetwork.Table_SR2_6.Co;
                        gExp = 0.5*s(4)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(4)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(4)*q*abs(q)/s(1)^3/s(2)^2,-s(4)*q*abs(q)/s(1)^2/s(2)^3,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        ZExp = [s(3)*(s(1)+s(2))/2/s(1)/s(2)];
                        dZdq = 0;
                        dZds = [-s(3)/2/s(1)^2,-s(3)/2/s(2)^2,(s(1)+s(2))/2/s(1)/s(2),0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR2_6 Pyramidal Diffuser, with wall';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR2_6};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Length(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR3_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H = s(1);
                        %             Wo = s(2);
                        %             W1 = s(3);
                        %             rho = s(4);
                        %             Ao = s(1)*s(2);
                        %             V = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR3_1.WoW1,DuctNetwork.Table_SR3_1.HW1};
                        Co_Table = DuctNetwork.Table_SR3_1.Co;
                        gExp = 0.5*s(4)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(4)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(4)*q*abs(q)/s(1)^3/s(2)^2,-s(4)*q*abs(q)/s(1)^2/s(2)^3,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        ZExp = [s(2)/s(3);s(1)/s(3)];
                        dZdq = [0;0];
                        dZds = [0,1/s(3),-s(2)/s(3)^2,0;1/s(3),0,-s(1)/s(3)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR3_1 Elbow, 90 Degree';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR3_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Width from Fan(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR4_1(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             Ho = s(1);
                        %             W = s(2);
                        %             H1 = s(3);
                        %             L = s(4);
                        %             rho = s(5);
                        %             Ao = s(1)*s(2);
                        %             A1 = s(3)*s(2);
                        %             Vo = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR4_1.tanTheta,DuctNetwork.Table_SR4_1.AoA1};
                        Co_Table = DuctNetwork.Table_SR4_1.Co;
                        gExp = 0.5*s(5)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(5)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(5)*q*abs(q)/s(1)^3/s(2)^2,-s(5)*q*abs(q)/s(1)^2/s(2)^3,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        if s(3) > s(1)
                            ZExp = [(s(3)-s(1))/2/s(4);s(1)/s(3)];
                            dZdq = [0;0];
                            dZds = [-1/2/s(4),0,1/2/s(4),-(s(3)-s(1))/2/s(4)^2,0;1/s(3),0,-s(1)/s(3)^2,0,0];
                        else
                            ZExp = [(s(1)-s(3))/2/s(4);s(1)/s(3)];
                            dZdq = [0;0];
                            dZds = [1/2/s(4),0,-1/2/s(4),-(s(1)-s(3))/2/s(4)^2,0;1/s(3),0,-s(1)/s(3)^2,0,0];
                        end
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR4_1 Transition, Rectangular, Symmetrical, Supply Air Systems';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR4_1};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Height(m)','Width(m)','Height from Fan(m)','Length(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = SR7_17(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        %             H1 = s(1);
                        %             W1 = s(2);
                        %             Ho = s(3);
                        %             Wo = s(4);
                        %             L = s(5);
                        %             rho = s(6);
                        %             Ao = s(3)*s(4);
                        %             A1 = s(1)*s(2);
                        %             V1 = q/s(1)/s(2);
                        GridVec = {DuctNetwork.Table_SR7_17.AoA1,DuctNetwork.Table_SR7_17.tanTheta};
                        Co_Table = DuctNetwork.Table_SR7_17.Co;
                        gExp = 0.5*s(6)*q*abs(q)/s(1)^2/s(2)^2;
                        dgdq = s(6)*abs(q)/s(1)^2/s(2)^2;
                        dgds = [-s(6)*q*abs(q)/s(1)^3/s(2)^2,-s(6)*q*abs(q)/s(1)^2/s(2)^3,0,0,0,0.5*q*abs(q)/s(1)^2/s(2)^2];
                        if abs(s(3)-s(1)) > abs(s(4)-s(2))% use abs(s(3)-s(1)) in Z expression
                            if s(3) > s(1)
                                ZExp = [s(3)*s(4)/s(1)/s(2);(s(3)-s(1))/2/s(5)];
                                dZdq = [0;0];
                                dZds = [-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0;-1/2/s(5),0,1/2/s(5),0,-(s(3)-s(1))/2/s(5)^2,0];
                            else
                                ZExp = [s(3)*s(4)/s(1)/s(2);(s(1)-s(3))/2/s(5)];
                                dZdq = [0;0];
                                dZds = [-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0;1/2/s(5),0,-1/2/s(5),0,-(s(1)-s(3))/2/s(5)^2,0];
                            end
                        else% use abs(s(4)-s(2)) in Z expression
                            if s(4) > s(2)
                                ZExp = [s(3)*s(4)/s(1)/s(2);(s(4)-s(2))/2/s(5)];
                                dZdq = [0;0];
                                dZds = [-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0;0,-1/2/s(5),0,1/2/s(5),-(s(4)-s(2))/2/s(5)^2,0];
                            else
                                ZExp = [s(3)*s(4)/s(1)/s(2);(s(2)-s(4))/2/s(5)];
                                dZdq = [0;0];
                                dZds = [-s(3)*s(4)/s(1)^2/s(2),-s(3)*s(4)/s(1)/s(2)^2,s(4)/s(1)/s(2),s(3)/s(1)/s(2),0,0;0,1/2/s(5),0,-1/2/s(5),-(s(2)-s(4))/2/s(5)^2,0];
                            end
                        end
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Co_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                        varargout{ii}=dP;
                    case 'dPdQ'
                        varargout{ii}=dPdQ;
                    case 'dPdS'
                        varargout{ii}=dPdS;
                    case 'Model_Description'
                        varargout{ii}='SR7_17 Pyramidal Diffuser at Centrifugal Fan Outlet';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.SR7_17};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:6};
                    case 'Parameter_Description'
                        varargout{ii}={'Height from Fan(m)','Width from Fan(m)','Height(m)','Width(m)','Length(m)','Density(kg/m^3)'};
                    case 'Parameter_LowerBound'
                        varargout{ii}=[1e-9,1e-9,1e-9,1e-9,1e-9,1e-6];
                    case 'Parameter_UpperBound'
                        varargout{ii}=[1e9,1e9,1e9,1e9,1e9,1e4];
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
    end
    
    properties (Constant)
        Table_CD3_6 = load('FittingData/CD3_6.mat');
        Table_CD3_9 = load('FittingData/CD3_9.mat');
        Table_CD3_17 = load('FittingData/CD3_17.mat');
        Table_CD6_1 = load('FittingData/CD6_1.mat');
        Table_CD9_1 = load('FittingData/CD9_1.mat');
        
        Table_CR3_1 = load('FittingData/CR3_1.mat');
        Table_CR3_3 = load('FittingData/CR3_3.mat');
        Table_CR3_6 = load('FittingData/CR3_6.mat');
        Table_CR3_17 = load('FittingData/CR3_17.mat');
        Table_CR6_1 = load('FittingData/CR6_1.mat');
        Table_CR6_4 = load('FittingData/CR6_4.mat');
        Table_CR9_1 = load('FittingData/CR9_1.mat');
        
        Table_ED1_1 = load('FittingData/ED1_1.mat');
        Table_ED1_3 = load('FittingData/ED1_3.mat');
        Table_ED5_1 = load('FittingData/ED5_1.mat');
        Table_ED5_2 = load('FittingData/ED5_2.mat');
        Table_ED5_3 = load('FittingData/ED5_3.mat');
        Table_ED5_4 = load('FittingData/ED5_4.mat');
        Table_ED7_2 = load('FittingData/ED7_2.mat');
        
        Table_ER4_3 = load('FittingData/ER4_3.mat');
        Table_ER5_1 = load('FittingData/ER5_1.mat');
        Table_ER5_3 = load('FittingData/ER5_3.mat');
        Table_ER5_4 = load('FittingData/ER5_4.mat');
        Table_ER5_5 = load('FittingData/ER5_5.mat');
        
        Table_SD5_1 = load('FittingData/SD5_1.mat');
        Table_SD5_3 = load('FittingData/SD5_3.mat');
        Table_SD5_9 = load('FittingData/SD5_9.mat');
        Table_SD5_18 = load('FittingData/SD5_18.mat');
        
        Table_SR2_1 = load('FittingData/SR2_1.mat');
        Table_SR2_3 = load('FittingData/SR2_3.mat');
        Table_SR2_5 = load('FittingData/SR2_5.mat');
        Table_SR2_6 = load('FittingData/SR2_6.mat');
        Table_SR3_1 = load('FittingData/SR3_1.mat');
        Table_SR4_1 = load('FittingData/SR4_1.mat');
        Table_SR5_1 = load('FittingData/SR5_1.mat');
        Table_SR5_13 = load('FittingData/SR5_13.mat');
        Table_SR5_15 = load('FittingData/SR5_15.mat');
        Table_SR5_14 = load('FittingData/SR5_14.mat');
        Table_SR7_17 = load('FittingData/SR7_17.mat');
        
        Table_NVT250 = load('FittingData/nvt250.mat');
    end
    
    methods (Static = true)
        function [dP, dPdQ, dPdS]=ED5_1(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ED5_1.QbQc,DuctNetwork.Table_ED5_1.AbAc,DuctNetwork.Table_ED5_1.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_1.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_ED5_1.QsQc,DuctNetwork.Table_ED5_1.AbAc,DuctNetwork.Table_ED5_1.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_1.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_1.QbQc,DuctNetwork.Table_ED5_1.AbAc,DuctNetwork.Table_ED5_1.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_1.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_1.QsQc,DuctNetwork.Table_ED5_1.AbAc,DuctNetwork.Table_ED5_1.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_1.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ED5_2(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ED5_2.QbQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_ED5_2.QsQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_2.QbQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_2.QsQc,DuctNetwork.Table_ED5_2.AbAc,DuctNetwork.Table_ED5_2.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ED5_2.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ED5_3(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        if s(3)<=0.25 %Dc <= 0.25m
                            Cb_Table = DuctNetwork.Table_ED5_3.Cb_part1;
                        else
                            Cb_Table = DuctNetwork.Table_ED5_3.Cb_part2;
                        end
                        GridVec = {DuctNetwork.Table_ED5_3.QbQc,DuctNetwork.Table_ED5_3.AbAc,DuctNetwork.Table_ED5_3.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cb_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        if s(3)<=0.25 %Dc <= 0.25m
                            Cs_Table = DuctNetwork.Table_ED5_3.Cs_part1;
                        else
                            Cs_Table = DuctNetwork.Table_ED5_3.Cs_part2;
                        end
                        GridVec = {DuctNetwork.Table_ED5_3.QsQc,DuctNetwork.Table_ED5_3.AbAc,DuctNetwork.Table_ED5_3.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cs_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        if s(3)<=0.25 %Dc <= 0.25m
                            Cb2_Table = DuctNetwork.Table_ED5_3.Cb2_part1;
                        else
                            Cb2_Table = DuctNetwork.Table_ED5_3.Cb2_part2;
                        end
                        GridVec = {DuctNetwork.Table_ED5_3.QbQc,DuctNetwork.Table_ED5_3.AbAc,DuctNetwork.Table_ED5_3.AsAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cb2_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        if s(3)<=0.25 %Dc <= 0.25m
                            Cs2_Table = DuctNetwork.Table_ED5_3.Cs2_part1;
                        else
                            Cs2_Table = DuctNetwork.Table_ED5_3.Cs2_part2;
                        end
                        GridVec = {DuctNetwork.Table_ED5_3.QsQc,DuctNetwork.Table_ED5_3.AbAc,DuctNetwork.Table_ED5_3.AsAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cs2_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ED5_4(q, s, Selection)
            % q = [Qb1,Qb2,Qc];
            % s = [Db1,Db2,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b1'
                        GridVec = {DuctNetwork.Table_ED5_4.QbQc,DuctNetwork.Table_ED5_4.AbAc,DuctNetwork.Table_ED5_4.AbAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ED5_4.Cb1, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 'b2'
                        GridVec = {DuctNetwork.Table_ED5_4.QbQc,DuctNetwork.Table_ED5_4.AbAc,DuctNetwork.Table_ED5_4.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ED5_4.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b1'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_4.QbQc,DuctNetwork.Table_ED5_4.AbAc,DuctNetwork.Table_ED5_4.AbAc};
                        ZExp = [q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ED5_4.Cb1_o, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 'b2'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_ED5_4.QbQc,DuctNetwork.Table_ED5_4.AbAc,DuctNetwork.Table_ED5_4.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ED5_4.Cb2_o, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ER5_1(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [H,Ws,Wb,Wc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(5)*(q(3)/(s(1)*s(4)))^2;
                dgdq =  [0,0,s(5)*q(3)/(s(1)*s(4))^2];
                dgds =  [-2/s(1),0,0,-2/s(4),1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ER5_1.QbQc,DuctNetwork.Table_ER5_1.AbAc,DuctNetwork.Table_ER5_1.AsAc};
                        ZExp = [q(2)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_1.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_ER5_1.QsQc,DuctNetwork.Table_ER5_1.AbAc,DuctNetwork.Table_ER5_1.AsAc};
                        ZExp = [q(1)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_1.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(5)*(q(2)/(s(1)*s(3)))^2;
                        dgdq =  [0,0,s(5)*q(2)/(s(1)*s(3))^2];
                        dgds =  [-2/s(1),0,-2/s(3),0,1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_ER5_1.QbQc,DuctNetwork.Table_ER5_1.AbAc,DuctNetwork.Table_ER5_1.AsAc};
                        ZExp = [q(2)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_1.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
                        dgdq =  [0,0,s(5)*q(1)/(s(1)*s(2))^2];
                        dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_ER5_1.QsQc,DuctNetwork.Table_ER5_1.AbAc,DuctNetwork.Table_ER5_1.AsAc};
                        ZExp = [q(1)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_1.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ER5_3(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Hs,Ws,Hb,Wb,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(5)*(q(3)/(s(1)*s(2)))^2;
                dgdq =  [0,0,s(5)*q(3)/(s(1)*s(2))^2];
                dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ER5_3.QbQc};
                        ZExp = [q(2)/q(3)];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2];
                        dZds = [0,0,0,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_ER5_3.QsQc};
                        ZExp = [q(1)/q(3)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2];
                        dZds = [0,0,0,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(5)*(q(2)/(s(3)*s(4)))^2;
                        dgdq =  [0,0,s(5)*q(2)/(s(3)*s(4))^2];
                        dgds =  [0,0,-2/s(3),-2/s(4),1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_ER5_3.QbQc};
                        ZExp = [q(2)/q(3)];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2];
                        dZds = [0,0,0,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
                        dgdq =  [0,0,s(5)*q(3)/(s(1)*s(2))^2];
                        dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_ER5_3.QsQc};
                        ZExp = [q(1)/q(3)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2];
                        dZds = [0,0,0,0,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=ER5_4(q, s, Selection)
            % q = [Qb1,Qb2,Qc];
            % s = [H,Wb1,Wb2,Wc,rho]
            gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
            dgdq =  [s(5)*q(1)/(s(1)*s(2))^2,0,0];
            dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
            switch Selection
                case 'b'
                    GridVec = {DuctNetwork.Table_ER5_4.AbAc};
                    ZExp = s(2)/s(4);
                    dZdq = [0,0,0];
                    dZds = [0,1/s(4),0,-s(2)/s(4)^2,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ER5_4.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=ER5_5(q, s, Selection)
            % q = [Qb1,Qb2,Qc];
            % s = [H,Wb1,Wb2,Wc,rho]
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(5)*(q(3)/(s(1)*s(4)))^2;
                dgdq =  [0,0,s(5)*q(3)/(s(1)*s(4))^2];
                dgds =  [-2/s(1),0,0,-2/s(4),1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ER5_5.QbQc,DuctNetwork.Table_ER5_5.AbAc};
                        ZExp = [q(1)/q(3);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ER5_5.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            else
                gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
                dgdq =  [0,0,s(5)*q(1)/(s(1)*s(2))^2];
                dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_ER5_5.QbQc,DuctNetwork.Table_ER5_5.AbAc};
                        ZExp = [q(1)/q(3);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_ER5_5.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SD5_1(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SD5_1.QbQc,DuctNetwork.Table_SD5_1.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_1.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_SD5_1.QsQc,DuctNetwork.Table_SD5_1.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_1.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_1.QbQc,DuctNetwork.Table_SD5_1.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_1.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_1.QsQc,DuctNetwork.Table_SD5_1.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_1.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SD5_3(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SD5_3.QbQc,DuctNetwork.Table_SD5_3.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_3.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_SD5_3.QsQc,DuctNetwork.Table_SD5_3.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_3.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_3.QbQc,DuctNetwork.Table_SD5_3.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_3.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_3.QsQc,DuctNetwork.Table_SD5_3.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_3.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SD5_9(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Ds,Db,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SD5_9.QbQc,DuctNetwork.Table_SD5_9.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_9.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_SD5_9.QsQc,DuctNetwork.Table_SD5_9.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_9.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(4)*(q(2)/(pi*s(2)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(2)/(pi*s(2)^2/4)^2];
                        dgds =  [0,-4/s(2),0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_9.QbQc,DuctNetwork.Table_SD5_9.AbAc};
                        ZExp = [q(2)/q(3);(s(2)/s(3))^2];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_9.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                        dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                        dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                        GridVec = {DuctNetwork.Table_SD5_9.QsQc,DuctNetwork.Table_SD5_9.AsAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_9.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SD5_18(q, s, Selection)
            % q = [Qb1,Qb2,Qc];
            % s = [Db1,Db2,Dc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
                dgdq =  [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
                dgds =  [0,0,-4/s(3),1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SD5_18.QbQc,DuctNetwork.Table_SD5_18.AbAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_18.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            else
                gExp =  0.5*s(4)*(q(1)/(pi*s(1)^2/4))^2;
                dgdq =  [0,0,s(4)*q(1)/(pi*s(1)^2/4)^2];
                dgds =  [-4/s(1),0,0,1/s(4)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SD5_18.QbQc,DuctNetwork.Table_SD5_18.AbAc};
                        ZExp = [q(1)/q(3);(s(1)/s(3))^2];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SD5_18.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SR5_1(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [H,Ws,Wb,Wc,rho];
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(5)*(q(3)/(s(1)*s(4)))^2;
                dgdq =  [0,0,s(5)*q(3)/(s(1)*s(4))^2];
                dgds =  [-2/s(1),0,0,-2/s(4),1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SR5_1.QbQc,DuctNetwork.Table_SR5_1.AbAc,DuctNetwork.Table_SR5_1.AsAc};
                        ZExp =[q(2)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq =[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds =[0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_1.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_SR5_1.QsQc,DuctNetwork.Table_SR5_1.AbAc,DuctNetwork.Table_SR5_1.AsAc};
                        ZExp = [q(1)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_1.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp =  0.5*s(5)*(q(2)/(s(1)*s(3)))^2;
                        dgdq =  [0,0,s(5)*q(2)/(s(1)*s(3))^2];
                        dgds =  [-2/s(1),0,-2/s(3),0,1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_SR5_1.QbQc,DuctNetwork.Table_SR5_1.AbAc,DuctNetwork.Table_SR5_1.AsAc};
                        ZExp =[q(2)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq =[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                        dZds =[0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_1.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
                        dgdq =  [0,0,s(5)*q(1)/(s(1)*s(2))^2];
                        dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                        GridVec = {DuctNetwork.Table_SR5_1.QsQc,DuctNetwork.Table_SR5_1.AbAc,DuctNetwork.Table_SR5_1.AsAc};
                        ZExp = [q(1)/q(3);s(3)/s(4);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                        dZds = [0,0,0,0,0;0,0,1/s(4),-s(3)/s(4)^2,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_1.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SR5_13(q, s, Selection)
            % q = [Qs,Qb,Qc];
            % s = [Hs,Ws,Hb,Wb,Hc,Wc,rho]
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp = 0.5*s(7)*(q(3)/(s(5)*s(6)))^2;
                dgdq = [0,0,s(7)*q(3)/(s(5)*s(6))^2];
                dgds = [0,0,0,0,-2/s(5),-2/s(6),1/s(7)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SR5_13.QbQc,DuctNetwork.Table_SR5_13.AbAc};
                        ZExp = [q(2)/q(3);(s(3)*s(4))/(s(5)*s(6))];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0,0,0;0,0,1/s(3),1/s(4),-1/s(5),-1/s(6),0]*(s(3)*s(4))/(s(5)*s(6));
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        GridVec = {DuctNetwork.Table_SR5_13.QsQc,DuctNetwork.Table_SR5_13.AsAc};
                        ZExp =[q(1)/q(3);(s(1)*s(2))/(s(5)*s(6))];
                        dZdq =[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds =[0,0,0,0,0,0,0;1/s(1),1/s(2),0,0,-1/s(5),-1/s(6),0]*(s(1)*s(2))/(s(5)*s(6));
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                end
            else
                switch Selection
                    case 'b'
                        gExp = 0.5*s(7)*(q(2)/(s(3)*s(4)))^2;
                        dgdq = [0,0,s(7)*q(2)/(s(3)*s(4))^2];
                        dgds = [0,0,-2/s(3),-2/s(4),0,0,1/s(7)]*gExp;
                        GridVec = {DuctNetwork.Table_SR5_13.QbQc,DuctNetwork.Table_SR5_13.AbAc};
                        ZExp = [q(2)/q(3);(s(3)*s(4))/(s(5)*s(6))];
                        dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0,0,0;0,0,1/s(3),1/s(4),-1/s(5),-1/s(6),0]*(s(3)*s(4))/(s(5)*s(6));
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    case 's'
                        gExp = 0.5*s(7)*(q(1)/(s(1)*s(2)))^2;
                        dgdq = [0,0,s(7)*q(1)/(s(1)*s(2))^2];
                        dgds = [-2/s(1),-2/s(2),0,0,0,0,1/s(7)]*gExp;
                        GridVec = {DuctNetwork.Table_SR5_13.QsQc,DuctNetwork.Table_SR5_13.AsAc};
                        ZExp =[q(1)/q(3);(s(1)*s(2))/(s(5)*s(6))];
                        dZdq =[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds =[0,0,0,0,0,0,0;1/s(1),1/s(2),0,0,-1/s(5),-1/s(6),0]*(s(1)*s(2))/(s(5)*s(6));
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS]=SR5_14(q, s, Selection)
            gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
            dgdq =  [s(5)*q(1)/(s(1)*s(2))^2,0,0];
            dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
            switch Selection
                case 'b'
                    GridVec = {DuctNetwork.Table_SR5_14.AbAc};
                    ZExp = s(2)/s(4);
                    dZdq = [0,0,0];
                    dZds = [0,1/s(4),0,-s(2)/s(4)^2,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SR5_14.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=SR5_15(q, s, Selection)
            % q = [Qb1,Qb2,Qc];
            % s = [H,Wb1,Wb2,Wc,rho]
            ModifiedTable = 1;
            if (ModifiedTable == 1)
                gExp =  0.5*s(5)*(q(3)/(s(1)*s(4)))^2;
                dgdq =  [0,0,s(5)*q(3)/(s(1)*s(4))^2];
                dgds =  [-2/s(1),0,0,-2/s(4),1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SR5_15.QbQc,DuctNetwork.Table_SR5_15.AbAc};
                        ZExp = [q(1)/q(3);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SR5_15.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            else
                gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
                dgdq =  [0,0,s(5)*q(1)/(s(1)*s(2))^2];
                dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
                switch Selection
                    case 'b'
                        GridVec = {DuctNetwork.Table_SR5_15.QbQc,DuctNetwork.Table_SR5_15.AbAc};
                        ZExp = [q(1)/q(3);s(2)/s(4)];
                        dZdq = [1/q(3),0,-q(1)/q(3)^2;0,0,0];
                        dZds = [0,0,0,0,0;0,1/s(4),0,-s(2)/s(4)^2,0];
                        [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.Table_SR5_15.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
                    otherwise
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
                end
            end
        end
        
        function [dP, dPdQ, dPdS] = Fan(q, s, c)
            v = s(1);
            if q<=0
                dP = -v.^2*exp(c(2));
                dPdQ = 0;
                dPdS = -2*v*exp(c(2));
            elseif q<v*exp(-c(1))
                m = c(3)/(log(q/v)+c(1));
                dP = -v^2*exp(m+c(2));
                dPdQ = -dP/c(3)/q*m^2;
                dPdS = dP/v*(m^2/c(3)+2);
            else %q>=v*exp(-c(1))
                dP = 0;
                dPdQ = 0;
                dPdS = 0;
            end
        end
        
        function [Eff, dEffdQ, dEffdS] = Fan_Efficiency(q, s, c, LnQ_V, Efficiency)
            X = log(q/s);
            Eff = interp1(LnQ_V,Efficiency,X,'linear','nearest');
            dEffdX = interp1(LnQ_V,[diff(Efficiency)./diff(LnQ_V);0],X,'previous',0);
            dEffdQ = dEffdX/q;
            dEffdS = -dEffdX/s;
        end
        
        function [f,dfdq,dfds] = Interp_Gradient(GridVector, InterpTable, Z, dZdq, dZds, g, dgdq, dgds)
            [Cf, dCfdZ] = DuctNetwork.InterpolationWithGradient(GridVector, InterpTable, Z);
            %             jac2 = DuctNetwork.Jacobian(@(x)DuctNetwork.InterpolationWithGradient(GridVector, InterpTable, x),Z,1e-6);
            %             jac3 = DuctNetwork.Jacobian(@(x)DuctNetwork.InterpolationWithGradient(GridVector, InterpTable, x),Z,-1e-6);
            %             if norm(jac2+jac3-2*dCfdZ)>1e-3
            %                 dCfdZ
            %                 jac
            %                 jac2
            %                 jac3
            %             end
            f = Cf*g;
            dfdq = g*dCfdZ*dZdq + Cf*dgdq;
            dfds = g*dCfdZ*dZds + Cf*dgds;
        end
        
        function [f,dfdq,dfds] = Double_Interp_Gradient(GridVector1, InterpTable1, GridVector2, InterpTable2, Z1, dZ1dq, dZ1ds,Z2, dZ2dq, dZ2ds, g, dgdq, dgds)
            [Cf1, dCf1dZ1] = DuctNetwork.InterpolationWithGradient(GridVector1, InterpTable1, Z1);
            [Cf2, dCf2dZ2] = DuctNetwork.InterpolationWithGradient(GridVector2, InterpTable2, Z2);
            f = g*Cf1*Cf2;
            dfdq = g*Cf2*dCf1dZ1*dZ1dq + g*Cf1*dCf2dZ2*dZ2dq + Cf1*Cf2*dgdq;
            dfds = g*Cf2*dCf1dZ1*dZ1ds + g*Cf1*dCf2dZ2*dZ2ds + Cf1*Cf2*dgds;
        end
        
        function [Cf, dCfdZ] = InterpolationWithGradient(GridVector, InterpTable, Z)
            range = @(x) max(x)-min(x);
            Z = reshape(Z,1,[]);
            NZ = length(Z);
            dCfdZ = zeros(1,NZ);
            n = cellfun(@length,GridVector);
            CfIdxInTable = cell(1,NZ);ZMesh = cell(1,NZ);IndexInMesh = cell(1,NZ);lambda = zeros(1,NZ);
            b_Interp = true(1,NZ);
            for ii = 1:NZ
                if Z(ii)>GridVector{ii}(end)
                    CfIdxInTable{ii}=n(ii)*[1;1;1];
                    ZMesh{ii}=GridVector{ii}(end)*[1;0;-1]+Z(ii)*[0;1;2];
                    b_Interp(ii) = false;
                    IndexInMesh{ii}=[1;2;3];
                elseif Z(ii)<GridVector{ii}(1)
                    CfIdxInTable{ii}=[1;1;1];
                    ZMesh{ii}=GridVector{ii}(1)*[-1;0;1]+Z(ii)*[2;1;0];
                    b_Interp(ii) = false;
                    IndexInMesh{ii}=[1;2;3];
                elseif any(Z(ii)==GridVector{ii})
                    ZLocation = find(Z(ii)==GridVector{ii},1);
                    if ZLocation==1
                        CfIdxInTable{ii}=[1;1;2];
                        ZMesh{ii}=reshape(GridVector{ii}(CfIdxInTable{ii}),[],1)-[1;0;0];
                        b_Interp(ii) = false;
                        IndexInMesh{ii}=[1;2;3];
                    elseif ZLocation==n(ii)
                        CfIdxInTable{ii}=n(ii)*[1;1;1]-[1;0;0];
                        ZMesh{ii}=reshape(GridVector{ii}(CfIdxInTable{ii}),[],1)+[0;0;1];
                        b_Interp(ii) = false;
                        IndexInMesh{ii}=[1;2;3];
                    else
                        CfIdxInTable{ii}=ZLocation*[1;1;1]+[-1;0;1];
                        ZMesh{ii}=reshape(GridVector{ii}(CfIdxInTable{ii}),[],1);
                        b_Interp(ii) = false;
                        IndexInMesh{ii}=[1;2;3];
                    end
                else
                    CfIdxInTable{ii}=find(Z(ii)<GridVector{ii},1)*[1;1]+[-1;0];
                    ZMesh{ii}=reshape(GridVector{ii}(CfIdxInTable{ii}),[],1);
                    lambda(ii) = (Z(ii)-ZMesh{ii}(1))/range(ZMesh{ii});
                    IndexInMesh{ii}=[1;2];
                end
            end
            CfMesh = InterpTable(CfIdxInTable{:});
            for ii=find(b_Interp)
                Index_1 = IndexInMesh; Index_1{ii}=1;
                Index_2 = IndexInMesh; Index_2{ii}=2;
                Index_3 = IndexInMesh; Index_3{ii}=3;
                CfMesh(Index_3{:}) = CfMesh(Index_2{:});
                CfMesh(Index_2{:}) = (1-lambda(ii))*CfMesh(Index_1{:}) + lambda(ii)*CfMesh(Index_3{:});
                ZMesh{ii} = [ZMesh{ii}(1); Z(ii);ZMesh{ii}(2)];
                IndexInMesh{ii}=[1;2;3];
            end
            Index_mid = num2cell(2*ones(1,NZ));
            Cf=CfMesh(Index_mid{:});
            for ii = 1:NZ
                Index_up = Index_mid; Index_up{ii}=3;
                Index_low = Index_mid; Index_low{ii}=1;
                dCfdZ(ii)=(CfMesh(Index_up{:})-Cf)/(ZMesh{ii}(3)-ZMesh{ii}(2))/2+(CfMesh(Index_low{:})-Cf)/(ZMesh{ii}(1)-ZMesh{ii}(2))/2;
            end
            
        end        
    end
end