classdef DuctNetworkSystem < matlab.System & matlab.system.mixin.Propagates ...
        & matlab.system.mixin.CustomIcon
    % Untitled Add summary here
    %
    % NOTE: When renaming the class name Untitled, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties, attributes,
    % and methods that you can implement for a System object in Simulink.
    
    % Public, tunable properties
    properties
    end
    
    properties
        n; % number of nodes without ground
        n_NodeDescription; % {'node1';'node2';} size n by 1, text description on each node
        P; %(Pa) array of size n by 1, pressure vector on each node
        
        b; % number of branches in the network
        b_FittingDescription; % {{'branch1_fit1';'branch1_fit2'};{'branch2_fit1'}}
        % cell array of size b by 1, text description on each branch
        b_Pdrop; % {{@dP_branch1_fit1;@dP_branch1_fit2};{@dP_branch2_fit1}}
        % cell array of size b by 1 for the pressure drop in terms of q and s
        % if there exist multiple Pdrop functions, use a cell array to store
        % each of them only record the original Pdrop function for the
        % fitting, handle the fitting direction in pressure drop calculation
        b_dPdQ; % {{@dPdQ_branch1_fit1;@dPdQ_branch1_fit2};{@dPdQ_branch2_fit1}}
        % cell array of size b by 1 for the partial pressure drop over q in terms of q and s
        % only record the original dPdQ function for the fitting, handle the
        % fitting direction in pressure drop calculation
        b_dPdS; % {{@dPdS_branch1_fit1;@dPdS_branch1_fit2};{@dPdS_branch2_fit1}}
        % cell array of size b by 1 for the partial pressure drop over s in terms of q and s
        % only record the original dPdS function for the fitting, handle the
        % fitting direction in pressure drop calculation
        b_Qidx; %{}
        % cell array of size b by 1 for the index of dependent Q of Pdrop, dPdQ and dPdS functions,
        % so use Q(abs(b_Qidx{b})).*sign(b_Qidx{b}) in these functions
        % Qidx{ii}[jj] can be negative, the sign represents the fitting
        % direction w.r.t. branch direction
        b_Sidx; % cell array of size b by 1 for the index of dependent S of Pdrop, dPdQ and dPdS functions,
        % so use S(b_Sidx{b}) in these functions
        Q; %(m^3/s) array of size b by 1, flow rate vector on each branch
        
        A; % Incidence matrix of size n by b
        % A(i,j)==1 means branch j leaves node i
        % A(i,j)==0 means branch j is not connected with node i
        % A(i,j)==-1 means branch j enters node i
        t; % dimension of null space of A (system's inner state)
        U; % null space basis of matrix A, of size b by t, AU=0, U'U=1
        X; % internal state of duct system, of size t by 1
        
        s; % number of parameters in the model
        s_ParamDescription; % cell array of size b by 1, text description on each parameter
        S; % array of size s by 1, parameter vector for whole system for identification
        
        d; % number of dampers
        d_Sidx; % array of parameters index for damper positions
        gdr; % number of terminals
        gdr_Bidx; % array of branch index for GDRs
        gdr_Nidx; % array of node index for GDRs
    end
    
    % Public, non-tunable properties
    properties(Nontunable)
        
    end
    
    properties(DiscreteState)
        
    end
    
    % Pre-computed constants
    properties(Access = private)
        
    end
    
    methods
        % Constructor
        function obj = DuctNetworkSystem(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods
        function Branch_Idx = AddBranch(obj, FromNode, ToNode, varargin)
            [~, b_Curr]=size(obj.A);
            Branch_Idx = b_Curr+1;
            SetNode(obj, FromNode, Branch_Idx, 1);
            SetNode(obj, ToNode, Branch_Idx, -1);
            [obj.n,obj.b]=size(obj.A);
            obj.U = 1e-1*null(obj.A);
            obj.t = size(obj.U,2);
            obj.X = ones(obj.t,1);
%             obj.X_lb = -1e12*ones(obj.t,1);
%             obj.X_ub = 1e12*ones(obj.t,1);
            obj.b_Pdrop{Branch_Idx,1}=cell(0,1);
            obj.b_dPdQ{Branch_Idx,1}=cell(0,1);
            obj.b_dPdS{Branch_Idx,1}=cell(0,1);
            obj.b_Qidx{Branch_Idx,1}=cell(0,1);
            obj.b_Sidx{Branch_Idx,1}=cell(0,1);
            obj.b_FittingDescription{Branch_Idx,1}=cell(0,1);
            if nargin>3
                obj.AddFitting(Branch_Idx,varargin{:});
            end
            function SetNode(obj, Node, Branch_Idx, a)
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
%                     obj.S_scale(Param_Idx(1:length(Param)),1) = Param;
%                     obj.S_ub(Param_Idx,1) = Model('Parameter_UpperBound');
%                     obj.S_lb(Param_Idx,1) = Model('Parameter_LowerBound');
%                     obj.s_m(Param_Idx,1) = Model('Is_Identified_Parameter');
%                     obj.s_MultiS(Param_Idx(1:length(Param)),1) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param,Param_Idx(1:length(Param)),'UniformOutput',false);
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
%                     obj.S_ub(Param_Idx,1) = Model('Parameter_UpperBound');
%                     obj.S_lb(Param_Idx,1) = Model('Parameter_LowerBound');
%                     obj.s_m(Param_Idx,1) = Model('Is_Identified_Parameter');
%                     obj.s_MultiS(Param_Idx(1:length(Param)),1) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param, Param_Idx(1:length(Param)),'UniformOutput',false);
                end
            elseif ischar(Model)
                obj.AddFitting(Branch,str2func(['Fittings.',Model]),Param);
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
%                 obj.S_scale(Param_Idx,1) = ParamValue;
%                 obj.S_ub(Param_Idx,1) = ParamUpperBound;
%                 obj.S_lb(Param_Idx,1) = ParamLowerBound;
%                 if nargin>=6
%                     obj.s_m(Param_Idx,1) = varargin{1};
%                     obj.s_MultiS{Param_Idx,1} = ParamValue*ones(1,varargin{1});
%                 else
%                     obj.s_m(Param_Idx,1) = 0;
%                     obj.s_MultiS{Param_Idx,1} = [];
%                 end
            end
        end
        
    end
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end
        
        function y = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            y = u;
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            
            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            
            % Set private and protected properties
            % obj.myproperty = s.myproperty;
            
            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        %% Simulink functions
        function ds = getDiscreteStateImpl(obj)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end
        
        function flag = isInputSizeLockedImpl(obj,index)
            % Return true if input size is not allowed to change while
            % system is running
            flag = true;
        end
        
        function out = getOutputSizeImpl(obj)
            % Return size for each output port
            out = [1 1];
            
            % Example: inherit size from first input port
            % out = propagatedInputSize(obj,1);
        end
        
        function icon = getIconImpl(obj)
            % Define icon for System block
            icon = mfilename('class'); % Use class name
            % icon = 'My System'; % Example: text icon
            % icon = {'My','System'}; % Example: multi-line text icon
            % icon = matlab.system.display.Icon('myicon.jpg'); % Example: image file icon
        end
    end
    
    methods(Static, Access = protected)
        %% Simulink customization functions
        function header = getHeaderImpl
            % Define header panel for System block dialog
            header = matlab.system.display.Header(mfilename('class'));
        end
        
        function group = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            group = matlab.system.display.Section(mfilename('class'));
        end
    end
end
