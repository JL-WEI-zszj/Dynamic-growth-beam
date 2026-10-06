function R = solve_generalized_kinematics_Ultimate_Vectorized_Video_dbug_C0_adaptive(user)

if nargin<1,user=struct();end
% ==================== 1. 用户基础设置 ====================
p=struct('h0',.01,'epsilon',.001,'nGrid',500,'Tmax',20, ...
    'beta',5e-8,'MassReg',1e-9,'RelTol',1e-7,'AbsTol',1e-9,'AuxAbsTol',1, ...
    'MaxStep',.05,'InitialStep',[],'MaxOrder',5,'MaxWallSeconds',10800, ...
    'PrintEverySeconds',15,'OutputDt',.02,'ExportTimes',[], ...
    'ShapeTimes',[],'PlotMode','Relative','MakeVideo',false,'PlotResults',true, ...
    'WriteFiles',true,'ExportShapeImage',true,'OutputDir','','FigureVisible','on','FPS',30,'PlaySpeed',1, ...
    'FigureSize',[900,600],'PlotBoxAspectRatio',1.5,'LegendFontSize',20, ...
    'ShapeFigureSize',[700,700],'ShapePlotBoxAspectRatio',1,'ShapeAxesMode','adaptive', ...
    'InitialAngle',0,'InitialVelocity',0,'F',[],'SelfTestOnly',false,'Stats','on');
assert(isstruct(user)&&isscalar(user),'输入必须为标量结构体。');
names=fieldnames(user);
for j=1:numel(names)
    if ~isfield(p,names{j}),error('BeamGeneral:Option','未知参数：%s',names{j});end
    p.(names{j})=user.(names{j});
end
if ~isfield(user,'ExportTimes'),p.ExportTimes=0:4:p.Tmax;end
if ~isfield(user,'ShapeTimes'),p.ShapeTimes=0:2:p.Tmax;end
p=validate_options(p);
if p.SelfTestOnly
    R=general_selftest(p);disp(R);assert(R.passed,'通用代数自检未通过。');
    fprintf('通用场Jacobian、消元、积分和边界自检通过；不等同于连续模型/网格验证。\n');
    return;
end
% p.F 可直接提供函数句柄；否则由下面原风格符号表达式自动求导。
if isempty(p.F),[F,expressions,constantCoefficients]=user_symbolic_fields();
else
    F=p.F;expressions=struct('source','user-supplied function handles');constantCoefficients=false;
end
validate_fields(F);
o=general_operators(p.nGrid);
f0=evaluate_fields(F,o,0);
[g0,v0]=initial_profiles(p,F,o);
[y0,yp0]=consistent_state(g0,v0,f0,p,o);
initialResidual=norm(general_residual(y0,yp0,f0,p,o),inf);
if ~isfinite(initialResidual)||initialResidual>1e-8*max(1,norm(yp0,inf))
    error('BeamGeneral:InitialResidual','初始动力学不相容：%g',initialResidual);
end
p.RunStamp=datestr(now,'HHMMSS');
if p.WriteFiles,p.OutputDir=prepare_output_directory(p);end
fprintf('\nh0=%g epsilon=%g nGrid=%d Tmax=%g；通用运动/生长，4N=%d\n', ...
    p.h0,p.epsilon,p.nGrid,p.Tmax,4*p.nGrid);
fprintf('beta=%g MassReg=%g；主变量容差 %g/%g，辅助量AbsTol=%g\n', ...
    p.beta,p.MassReg,p.RelTol,p.AbsTol,p.AuxAbsTol);
[t,Y,meta]=integrate_general(p,F,o,y0,yp0,constantCoefficients);
k=general_postprocess(t,Y,p,F,o);
R=struct('t',t,'S',o.S,'S_angle',o.Sc,'g',Y(1:o.n,:)', ...
    'g_tau',Y(o.n+1:2*o.n,:)','g_tautau',k.a,'g_S',k.gS, ...
    'f_vals',k.f,'lambda0',k.lambda,'x0_vals',k.x,'z0_vals',k.z, ...
    'completed',meta.completed,'reason',meta.reason,'exception',meta.exception, ...
    'wall_solve_seconds',meta.wall_solve_seconds,'integrator',meta, ...
    'parameters',p,'fields',F,'field_expressions',expressions, ...
    'algebraic_residual',k.residual,'max_abs_strain',max(abs(k.f(:)./k.lambda(:)-1)), ...
    'initial_residual',initialResidual,'output_kind','requested_time_samples', ...
    'spatial_scheme','cell-centred Green kernel and ghost-point differences');
if p.WriteFiles,write_paper_outputs(R,p);end
if p.PlotResults||p.MakeVideo,plot_paper_outputs(R,p);end
fprintf('状态=%s；实际结束tau=%.12g；求解墙钟=%.2f秒；输出点=%d\n', ...
    R.reason,R.t(end),R.wall_solve_seconds,numel(R.t));
fprintf('采样点 max|f/lambda0-1|=%.6g；不能据完成状态判断网格收敛。\n',R.max_abs_strain);
if p.WriteFiles,fprintf('论文绘图数据目录：%s\n',p.OutputDir);end
end

function [F,expressions,constantCoefficients]=user_symbolic_fields()
% ==================== 2. 原风格：修改这五个表达式 ====================
% S、t均为无量纲变量；cx/cz只能依赖t。需要Symbolic Math Toolbox。
% 使用数值函数句柄 p.F 时完全不需要符号工具箱。
syms t_sym S_sym real
cx_sym = 0.5 * (cos(2*sym(pi)*0.1 * (t_sym - (1-exp(-2*t_sym))/2)) - 1);
cz_sym = 0.5 * sin(2*sym(pi)*0.1 * (t_sym - (1-exp(-2*t_sym))/2));
Nl10_sym = 1 + 0 * S_sym + 0 * t_sym;
Nl11_sym = 0 * S_sym + 0 * t_sym;
Nl12_sym = 0 * S_sym + 0 * t_sym;
% 示例：Nl10_sym=1+0.02*S_sym*sin(t_sym);
%       Nl11_sym=0.01*(1+S_sym)*(1-cos(t_sym));
% 若初始Nl11(1,0)或其时间导数非零，应设置相容的InitialAngle/Velocity。
if has(cx_sym,S_sym)||has(cz_sym,S_sym)
    error('BeamGeneral:BaseMotion','基座位移cx、cz不能依赖S。');
end
F.cx=matlabFunction(cx_sym,'Vars',{t_sym});
F.cz=matlabFunction(cz_sym,'Vars',{t_sym});
F.cxddot=matlabFunction(diff(cx_sym,t_sym,2),'Vars',{t_sym});
F.czddot=matlabFunction(diff(cz_sym,t_sym,2),'Vars',{t_sym});
F.Nl10=matlabFunction(Nl10_sym,'Vars',{S_sym,t_sym});
F.Nl10_S=matlabFunction(diff(Nl10_sym,S_sym),'Vars',{S_sym,t_sym});
F.Nl10_t=matlabFunction(diff(Nl10_sym,t_sym),'Vars',{S_sym,t_sym});
F.Nl10_tt=matlabFunction(diff(Nl10_sym,t_sym,2),'Vars',{S_sym,t_sym});
F.Nl11=matlabFunction(Nl11_sym,'Vars',{S_sym,t_sym});
F.Nl11_S=matlabFunction(diff(Nl11_sym,S_sym),'Vars',{S_sym,t_sym});
F.Nl11_t=matlabFunction(diff(Nl11_sym,t_sym),'Vars',{S_sym,t_sym});
F.Nl12=matlabFunction(Nl12_sym,'Vars',{S_sym,t_sym});
expressions=struct('cx',char(cx_sym),'cz',char(cz_sym),'Nl10',char(Nl10_sym), ...
    'Nl11',char(Nl11_sym),'Nl12',char(Nl12_sym));
% 只有符号表达式能证明系数不随时间变化时才缓存，任意数值F不作猜测。
constantCoefficients=~any(has([Nl10_sym,Nl11_sym,Nl12_sym, ...
    diff(cx_sym,t_sym,2),diff(cz_sym,t_sym,2)],t_sym));
end

function p=validate_options(p)
for name={'h0','epsilon','Tmax','RelTol','AbsTol','AuxAbsTol','MaxStep', ...
        'MaxWallSeconds','PrintEverySeconds','OutputDt','FPS','PlaySpeed'}
    validateattributes(p.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'beta','MassReg'}
    validateattributes(p.(name{1}),{'numeric'},{'scalar','real','finite','nonnegative'});
end
for name={'MakeVideo','PlotResults','WriteFiles','ExportShapeImage','SelfTestOnly'}
    validateattributes(p.(name{1}),{'numeric','logical'},{'scalar','real','finite','binary'});
end
validateattributes(p.nGrid,{'numeric'},{'scalar','integer','>=',6,'finite'});
validateattributes(p.MaxOrder,{'numeric'},{'scalar','integer','>=',1,'<=',5});
if ~isempty(p.InitialStep)
    validateattributes(p.InitialStep,{'numeric'},{'scalar','real','finite','positive'});
end
for name={'ExportTimes','ShapeTimes'}
    v=p.(name{1});
    validateattributes(v,{'numeric'},{'real','finite','nonnegative','<=',p.Tmax});
    if ~isempty(v),validateattributes(v,{'numeric'},{'vector'});end
    p.(name{1})=sort(unique(v(:)'));
end
p.PlotMode=validatestring(p.PlotMode,{'Relative','Absolute'});
p.FigureVisible=validatestring(p.FigureVisible,{'on','off'});
p.ShapeAxesMode=validatestring(p.ShapeAxesMode,{'adaptive','fixed'});
p.Stats=validatestring(p.Stats,{'on','off'});
for name={'FigureSize','ShapeFigureSize'}
    wh=p.(name{1});
    validateattributes(wh,{'numeric'}, ...
        {'vector','numel',2,'real','finite','integer','positive'});
    p.(name{1})=double(wh(:)');
    if any(p.(name{1})<[640,480])
        error('BeamGeneral:FigureSize','%s 至少为 [640,480] 像素，以容纳坐标文字。',name{1});
    end
end
validateattributes(p.PlotBoxAspectRatio,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.ShapePlotBoxAspectRatio,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(p.LegendFontSize,{'numeric'},{'scalar','real','finite','positive'});
% 旧设置中的矩形窗口/比例也统一为正方形；保留选项名以兼容原调用。
p.ShapeFigureSize=repmat(max(p.ShapeFigureSize),1,2);
p.ShapePlotBoxAspectRatio=1;
assert(ischar(p.OutputDir)||(isstring(p.OutputDir)&&isscalar(p.OutputDir)), ...
    'OutputDir必须为路径字符串。');p.OutputDir=char(p.OutputDir);
if p.MakeVideo&&(~p.WriteFiles||strcmp(p.FigureVisible,'off'))
    error('BeamGeneral:VideoOptions','视频需要 WriteFiles=true、FigureVisible=on。');
end
end

function validate_fields(F)
names={'cx','cz','cxddot','czddot','Nl10','Nl10_S','Nl10_t','Nl10_tt', ...
    'Nl11','Nl11_S','Nl11_t','Nl12'};
if ~isstruct(F)||~isscalar(F)||~all(isfield(F,names))
    error('BeamGeneral:Fields','F必须包含cx/cz及其二阶导数、Nl10及S/t/tt导数、Nl11及S/t导数、Nl12。');
end
for j=1:numel(names)
    if ~isa(F.(names{j}),'function_handle')
        error('BeamGeneral:Fields','F.%s必须是函数句柄。',names{j});
    end
end
end

function o=general_operators(n)
ds=1/n;Sc=((1:n)'-.5)*ds;
D2=spdiags(repmat([1,-2,1],n,1),-1:1,n,n)/ds^2;
D2(1,1)=-3/ds^2;D2(n,n)=-1/ds^2;
D1=spdiags(repmat([-1,0,1],n,1),-1:1,n,n)/(2*ds);
D1(1,1)=1/(2*ds);D1(n,n)=1/(2*ds);
Lp=spdiags(repmat([-1,2,-1],n,1),-1:1,n,n);Lp(1,1)=1;Lp(n,n)=3;
% H=ds*(1-max(Sc,Sc'))=ds^2*inv(Lp)，积分边界与角度边界不同。
o=struct('n',n,'ds',ds,'Sc',Sc,'S',(0:n)'*ds,'D1',D1,'D2',D2, ...
    'Lp',Lp,'I',speye(n),'Z',sparse(n,n));
end

function f=evaluate_fields(F,o,t)
S=o.Sc;
f.l=field_vector(F.Nl10(S,t),o.n,'Nl10');
f.ls=field_vector(F.Nl10_S(S,t),o.n,'Nl10_S');
f.lt=field_vector(F.Nl10_t(S,t),o.n,'Nl10_t');
f.ltt=field_vector(F.Nl10_tt(S,t),o.n,'Nl10_tt');
f.m=field_vector(F.Nl11(S,t),o.n,'Nl11');
f.ms=field_vector(F.Nl11_S(S,t),o.n,'Nl11_S');
f.m2=field_vector(F.Nl12(S,t),o.n,'Nl12');
% 边界值必须在S=1上求值，不能拿最后一个单元中心代替。
f.B=-field_vector(F.Nl11(1,t),1,'Nl11(1,t)');
f.Bt=-field_vector(F.Nl11_t(1,t),1,'Nl11_t(1,t)');
endStretch=field_vector(F.Nl10([0;1],t),2,'Nl10(endpoints,t)');
if any(f.l<=0)||any(endStretch<=0)
    error('BeamGeneral:Stretch','Nl10必须保持正数；tau=%.16g。',t);
end
f.bx=(1-S)*field_vector(F.cxddot(t),1,'cxddot');
f.bz=(1-S)*field_vector(F.czddot(t),1,'czddot');
end

function v=field_vector(v,n,name)
if ~isnumeric(v)||~isreal(v)||any(~isfinite(v(:)))||(~isscalar(v)&&numel(v)~=n)
    error('BeamGeneral:FieldValue','%s必须返回有限实数标量或与S等长的向量。',name);
end
if isscalar(v),v=repmat(double(v),n,1);else,v=double(v(:));end
end

function [g,v]=initial_profiles(p,F,o)
values={p.InitialAngle,p.InitialVelocity};boundary=[-F.Nl11(1,0),-F.Nl11_t(1,0)];
result=cell(1,2);
for j=1:2
    value=values{j};
    if isa(value,'function_handle')
        result{j}=field_vector(value(o.Sc),o.n,'initial profile');
        root=field_vector(value(0),1,'initial profile at root');
        % 五点单边差分仅检查输入函数在自由端的边界相容性。
        d=1e-4;tipSamples=field_vector(value(1-(0:4)'*d),5,'initial profile near tip');
        tipDerivative=[25,-48,36,-16,3]*tipSamples/(12*d);
        if abs(root)>1e-8||abs(tipDerivative-boundary(j))>1e-5*max(1,abs(boundary(j)))
            error('BeamGeneral:InitialBoundary','初始函数必须满足根部值0及自由端导数-Nl11或-Nl11_t。');
        end
    elseif isnumeric(value)&&isscalar(value)
        if value~=0||abs(boundary(j))>1e-12
            error('BeamGeneral:InitialBoundary','零初值与当前端部生长不相容；请提供相容的InitialAngle/InitialVelocity函数。');
        end
        result{j}=zeros(o.n,1);
    else
        result{j}=field_vector(value,o.n,'initial profile vector');
        % 向量为单元中心值；根部/自由端条件由虚节点施加。
        % 无法由有限个内部样本严格确认连续初始场的边界相容性。
        warning('BeamGeneral:InitialSamples','初值向量按单元中心解释；请自行确认其连续场满足边界条件。');
    end
end
g=result{1};v=result{2};
end

function [uS,uSS]=angle_derivatives(u,B,o)
% g_left=-g_1, g_right=g_N+ds*B。B=-Nl11(1,t)。
uS=o.D1*u;uSS=o.D2*u;
uS(end)=uS(end)+B/2;uSS(end)=uSS(end)+B/o.ds;
end

function b=bending_term(g,v,f,p,o)
[uS,uSS]=angle_derivatives(g+p.beta*v,f.B+p.beta*f.Bt,o);
b=(8*p.h0^2/(3*p.epsilon))./f.l.^4.* ...
    (-2*(f.m+uS).*f.ls+f.l.*(f.ms+uSS));
end

function K=bending_jacobian(f,p,o)
K=spdiags((8*p.h0^2/(3*p.epsilon))./f.l.^4,0,o.n,o.n)* ...
    (-2*spdiags(f.ls,0,o.n,o.n)*o.D1+spdiags(f.l,0,o.n,o.n)*o.D2);
end

function r=general_residual(y,yp,f,p,o)
n=o.n;g=y(1:n);v=y(n+1:2*n);px=y(2*n+1:3*n);pz=y(3*n+1:4*n);
a=yp(n+1:2*n);c=cos(g);s=sin(g);
qx=c.*(f.ltt-f.l.*v.^2)-s.*(2*f.lt.*v+f.l.*a);
qz=s.*(f.ltt-f.l.*v.^2)+c.*(2*f.lt.*v+f.l.*a);
r=[yp(1:n)-v;-bending_term(g,v,f,p,o)-s.*(px+f.bx)+c.*(pz+f.bz)+ ...
    (p.MassReg/p.epsilon)*a;o.Lp*px-o.ds^2*qx;o.Lp*pz-o.ds^2*qz];
end

function [Jy,Jyp]=general_jacobian(y,yp,f,p,o)
n=o.n;g=y(1:n);v=y(n+1:2*n);px=y(2*n+1:3*n);pz=y(3*n+1:4*n);
a=yp(n+1:2*n);c=cos(g);s=sin(g);d=o.ds^2;
qx=c.*(f.ltt-f.l.*v.^2)-s.*(2*f.lt.*v+f.l.*a);
qz=s.*(f.ltt-f.l.*v.^2)+c.*(2*f.lt.*v+f.l.*a);
C=spdiags(c,0,n,n);S=spdiags(s,0,n,n);I=o.I;Z=o.Z;K=bending_jacobian(f,p,o);
Jy=[Z,-I,Z,Z; ...
    -K+spdiags(-c.*(px+f.bx)-s.*(pz+f.bz),0,n,n),-p.beta*K,-S,C; ...
    d*spdiags(qz,0,n,n),2*d*spdiags(c.*f.l.*v+s.*f.lt,0,n,n),o.Lp,Z; ...
    -d*spdiags(qx,0,n,n),2*d*spdiags(s.*f.l.*v-c.*f.lt,0,n,n),Z,o.Lp];
Jyp=[I,Z,Z,Z;Z,(p.MassReg/p.epsilon)*I,Z,Z; ...
    Z,d*spdiags(f.l.*s,0,n,n),Z,Z;Z,-d*spdiags(f.l.*c,0,n,n),Z,Z];
end

function [y,yp]=consistent_state(g,v,f,p,o)
n=o.n;c=cos(g);s=sin(g);C=spdiags(c,0,n,n);S=spdiags(s,0,n,n);d=o.ds^2;
qxv=c.*(f.ltt-f.l.*v.^2)-s.*(2*f.lt.*v);
qzv=s.*(f.ltt-f.l.*v.^2)+c.*(2*f.lt.*v);
A=[(p.MassReg/p.epsilon)*o.I,-S,C; ...
    d*spdiags(f.l.*s,0,n,n),o.Lp,o.Z; ...
    -d*spdiags(f.l.*c,0,n,n),o.Z,o.Lp];
u=A\[bending_term(g,v,f,p,o)+s.*f.bx-c.*f.bz;d*qxv;d*qzv];
if any(~isfinite(u))
    error('BeamGeneral:NonfiniteState','稀疏消元产生非有限数值，请检查参数、场函数和网格。');
end
y=[g;v;u(n+1:end)];yp=[v;u(1:n);zeros(2*n,1)];
end

function [t,Y,meta]=integrate_general(p,F,o,y0,yp0,constantCoefficients)
requested=unique([0:p.OutputDt:p.Tmax,p.ExportTimes,p.ShapeTimes,p.Tmax]);
if numel(requested)==2,requested=unique([requested,p.Tmax/2]);end
sampleT=zeros(numel(requested)+1,1);sampleY=zeros(2*o.n,numel(requested)+1);
count=1;sampleY(:,1)=y0(1:2*o.n);lastT=0;trialT=0;lastPrint=0;
nres=0;njac=0;nout=0;reason='';exception='';cachedT=NaN;cachedF=[];
clock0=tic;
opts=odeset('RelTol',p.RelTol,'AbsTol',[p.AbsTol*ones(2*o.n,1);p.AuxAbsTol*ones(2*o.n,1)], ...
    'MaxStep',p.MaxStep,'MaxOrder',p.MaxOrder,'Jacobian',@jacfun, ...
    'OutputFcn',@monitor,'Refine',1,'Stats',p.Stats);
if ~isempty(p.InitialStep),opts=odeset(opts,'InitialStep',p.InitialStep);end
try
    [t,rows]=ode15i(@resfun,requested,y0,yp0,opts);Y=rows(:,1:2*o.n)';
catch ME
    exception=ME.message;if isempty(reason),reason=['solver_error: ',ME.identifier];end
    t=sampleT(1:count);Y=sampleY(:,1:count);
    warning('BeamGeneral:Stopped','%s；仅保留至tau=%g的输出。',ME.message,t(end));
end
completed=t(end)>=p.Tmax-1e-10*max(1,p.Tmax);
if completed,reason='completed';elseif isempty(reason),reason='solver_returned_before_Tmax';end
meta=struct('completed',completed,'reason',reason,'exception',exception, ...
    'wall_solve_seconds',toc(clock0),'residual_calls',nres,'jacobian_calls',njac, ...
    'output_callbacks',nout,'last_trial_time',trialT,'method','ode15i-general-sparse');
    function f=coefficients(tt)
        if isempty(cachedF)||(~constantCoefficients&&~isequal(tt,cachedT))
            cachedF=evaluate_fields(F,o,tt);cachedT=tt;
        end
        f=cachedF;
    end
    function r=resfun(tt,y,yp)
        nres=nres+1;trialT=tt;check_budget();r=general_residual(y,yp,coefficients(tt),p,o);
    end
    function [A,B]=jacfun(tt,y,yp)
        njac=njac+1;check_budget();[A,B]=general_jacobian(y,yp,coefficients(tt),p,o);
    end
    function check_budget()
        elapsed=toc(clock0);
        if elapsed-lastPrint>=p.PrintEverySeconds
            fprintf('wall=%.1fs | tau saved=%.9g | trial=%.9g\n',elapsed,lastT,trialT);lastPrint=elapsed;
        end
        if elapsed>=p.MaxWallSeconds
            reason='wall_time_limit';error('BeamGeneral:WallLimit','达到MaxWallSeconds限时。');
        end
    end
    function status=monitor(tt,yy,flag)
        status=0;
        if isempty(flag)&&~isempty(tt)
            nout=nout+1;
            for k=1:numel(tt)
                if tt(k)>sampleT(count)
                    count=count+1;sampleT(count)=tt(k);sampleY(:,count)=yy(1:2*o.n,k);
                end
            end
            lastT=tt(end);
            if toc(clock0)>=p.MaxWallSeconds,reason='wall_time_limit';status=1;end
        end
    end
end

function k=general_postprocess(t,Y,p,F,o)
nt=numel(t);n=o.n;
k.a=zeros(nt,n);k.gS=k.a;k.f=k.a;k.lambda=k.a;k.x=zeros(nt,n+1);k.z=k.x;k.residual=zeros(nt,1);
for j=1:nt
    g=Y(1:n,j);v=Y(n+1:2*n,j);f=evaluate_fields(F,o,t(j));
    % 与推进时相同的方程、beta、边界和MassReg恢复加速度。
    [y,yp]=consistent_state(g,v,f,p,o);px=y(2*n+1:3*n);pz=y(3*n+1:4*n);
    c=cos(g);s=sin(g);[gS,~]=angle_derivatives(g,f.B,o);
    Bterm=-2*p.h0*(f.m-3*gS).*(f.m+gS)+f.l.*(3*f.m+2*p.h0*f.m2+3*gS);
    % 原式(31)的代数等价形式，避免先除epsilon再乘epsilon。
    stretch=f.l-(p.epsilon/8)*f.l.^2.*(c.*(px+f.bx)+s.*(pz+f.bz)) ...
        +(p.h0/3)*Bterm./f.l;
    if any(~isfinite(stretch))
        error('BeamGeneral:NonfiniteStretch','后处理伸长非有限，tau=%.16g。',t(j));
    end
    k.a(j,:)=yp(n+1:2*n)';k.gS(j,:)=gS';k.f(j,:)=stretch';k.lambda(j,:)=f.l';
    k.x(j,2:end)=o.ds*cumsum(stretch.*c)';k.z(j,2:end)=o.ds*cumsum(stretch.*s)';
    if strcmp(p.PlotMode,'Absolute')
        k.x(j,:)=k.x(j,:)+field_vector(F.cx(t(j)),1,'cx');
        k.z(j,:)=k.z(j,:)+field_vector(F.cz(t(j)),1,'cz');
    end
    k.residual(j)=norm(general_residual(y,yp,f,p,o),inf);
end
end

function outdir=prepare_output_directory(p)
outdir=p.OutputDir;
if isempty(outdir)
    base=fullfile(pwd,sprintf('beam_general_h%g_e%g_N%d_%s',p.h0,p.epsilon,p.nGrid,datestr(now,'yyyymmdd_HHMMSS')));
    outdir=base;j=1;
    while exist(outdir,'file')||exist(outdir,'dir'),outdir=sprintf('%s_%03d',base,j);j=j+1;end
end
if ~exist(outdir,'dir'),mkdir(outdir);end
[~,info]=fileattrib(outdir);outdir=info.Name;
planned={sprintf('Rod_TipZ_Displacement_%s.csv',p.RunStamp)};
for tt=p.ExportTimes,planned{end+1}=sprintf('Rod_Shape_tau=%.15g_%s.txt',tt,p.RunStamp);end %#ok<AGROW>
if p.MakeVideo,planned{end+1}=sprintf('Rod_Evolution_%s.mp4',p.RunStamp);end
for j=1:numel(planned)
    if exist(fullfile(outdir,planned{j}),'file')
        error('BeamGeneral:ExistingOutput','拒绝覆盖已有文件%s；请换OutputDir。',planned{j});
    end
end
end

function T=general_selftest(p)
old=rng;cleanup=onCleanup(@()rng(old));rng(197); %#ok<NASGU>
o=general_operators(20);n=o.n;ds=o.ds;F=selftest_fields();f=evaluate_fields(F,o,.37);
y=.1*randn(4*n,1);yp=.2*randn(4*n,1);u=randn(4*n,1);w=randn(4*n,1);d=1e-30;
[J,B]=general_jacobian(y,yp,f,p,o);
T.jac_y_error=norm(imag(general_residual(y+1i*d*u,yp,f,p,o))/d-J*u)/max(1,norm(J*u));
T.jac_yp_error=norm(imag(general_residual(y,yp+1i*d*w,f,p,o))/d-B*w)/max(1,norm(B*w));
H=ds*(1-max(o.Sc,o.Sc'));
T.integral_inverse_error=norm(H-ds^2*full(o.Lp\o.I),inf)/norm(H,inf);
T.constant_integral_error=max(abs(H*ones(n,1)-(1-o.Sc.^2)/2));
g=y(1:n);v=y(n+1:2*n);c=cos(g);s=sin(g);a=yp(n+1:2*n);
qxv=c.*(f.ltt-f.l.*v.^2)-s.*(2*f.lt.*v);
qzv=s.*(f.ltt-f.l.*v.^2)+c.*(2*f.lt.*v);
y(2*n+1:3*n)=H*(qxv-f.l.*s.*a);y(3*n+1:end)=H*(qzv+f.l.*c.*a);
r=general_residual(y,yp,f,p,o);
M=(H.*cos(g-g')).*f.l'+(p.MassReg/p.epsilon)*eye(n);
dense=M*a-bending_term(g,v,f,p,o)-s.*(H*qxv+f.bx)+c.*(H*qzv+f.bz);
T.elimination_error=norm(r(n+1:2*n)-dense)/max(1,norm(dense));
[yy,pp]=consistent_state(g,v,f,p,o);
T.initial_residual=norm(general_residual(yy,pp,f,p,o),inf);
[d1,d2]=angle_derivatives(.23*o.Sc,.23,o);
T.nonzero_boundary_linear_error=max([norm(d1-.23,inf),norm(d2,inf)]);
% 常生长诱导曲率：g=-m*S、Nl10=1、Nl11=m、无载，应为静态解。
Q=F;Q.Nl10=@(S,t)1;Q.Nl10_S=@(S,t)0;Q.Nl10_t=@(S,t)0;Q.Nl10_tt=@(S,t)0;
Q.Nl11=@(S,t).1;Q.Nl11_S=@(S,t)0;Q.Nl11_t=@(S,t)0;Q.Nl12=@(S,t)0;
Q.cxddot=@(t)0;Q.czddot=@(t)0;
fq=evaluate_fields(Q,o,0);[yq,pq]=consistent_state(-.1*o.Sc,zeros(n,1),fq,p,o);
T.static_curvature_acceleration=norm(pq(n+1:2*n),inf);
k=general_postprocess(0,yq(1:2*n),p,Q,o);T.static_curvature_stretch_error=norm(k.f-1,inf);
T.passed=max([T.jac_y_error,T.jac_yp_error,T.integral_inverse_error,T.elimination_error, ...
    T.initial_residual,T.nonzero_boundary_linear_error,T.static_curvature_acceleration, ...
    T.static_curvature_stretch_error])<1e-8&&abs(T.constant_integral_error-ds^2/8)<1e-12;
end

function F=selftest_fields()
F.cx=@(t).03*sin(t);F.cz=@(t).04*(1-cos(t));
F.cxddot=@(t)-.03*sin(t);F.czddot=@(t).04*cos(t);
F.Nl10=@(S,t)1+.1*S+.03*sin(t);F.Nl10_S=@(S,t).1;
F.Nl10_t=@(S,t).03*cos(t);F.Nl10_tt=@(S,t)-.03*sin(t);
F.Nl11=@(S,t).02*(1+S).*sin(t);F.Nl11_S=@(S,t).02*sin(t);
F.Nl11_t=@(S,t).02*(1+S).*cos(t);F.Nl12=@(S,t).01*cos(S+t);
end

%% ==================== 原版论文数据输出格式 ====================
function write_paper_outputs(R,p)
% CSV/TXT 的文件名、列顺序、表头和分隔符沿用原程序。
% ExportTimes 已加入求解器输出时刻；直接取该时刻的解，不作外推。
if ~p.WriteFiles, return; end
if ~exist(p.OutputDir,'dir'), mkdir(p.OutputDir); end
[exportTimes,exportRows] = paper_sample_rows(R.t,p.ExportTimes);
csvName = fullfile(p.OutputDir,sprintf('Rod_TipZ_Displacement_%s.csv',p.RunStamp));
txtNames = cell(numel(exportTimes),1);
for k = 1:numel(exportTimes)
    txtNames{k} = fullfile(p.OutputDir,sprintf('Rod_Shape_tau=%.15g_%s.txt', ...
        exportTimes(k),p.RunStamp));
end
% 不覆盖以前计算的数据；同一秒重复运行时可换 OutputDir。
paper_assert_new_files([{csvName};txtNames]);
fid = fopen(csvName,'w');
if fid < 0, error('beam:OutputOpen','无法创建输出文件：%s',csvName); end
closeFile = onCleanup(@() fclose(fid));
if strcmpi(p.PlotMode,'Relative')
    fprintf(fid,'tau,Tip_Z_Dimensionless_Relative\n');
else
    fprintf(fid,'tau,Tip_Z_Dimensionless_Absolute\n');
end
clear closeFile;
writematrix([R.t(:),R.z0_vals(:,end)],csvName,'WriteMode','append');
fprintf('末端 Z 位移数据已保存至: %s\n',csvName);
for k = 1:numel(exportTimes)
    row = exportRows(k);
    writematrix([R.x0_vals(row,:).',R.z0_vals(row,:).'],txtNames{k}, ...
        'Delimiter','tab');
    fprintf('  -> 已导出 tau = %.15g 的两列构型数据：%s\n', ...
        exportTimes(k),txtNames{k});
end
if numel(exportTimes) < numel(p.ExportTimes)
    fprintf('仅导出实际到达区间 [%.15g, %.15g] 内的构型，未生成后续时刻数据。\n', ...
        R.t(1),R.t(end));
end
end

function plot_paper_outputs(R,p)
% 构型图窗口/图框固定为 1:1，坐标范围自适应；末端位移图为 1.5:1。
% 横纵坐标始终等单位长度，不通过竖向拉伸构型来消除留白。
if ~p.PlotResults && ~p.MakeVideo, return; end
if p.MakeVideo && ~p.WriteFiles
    error('beam:VideoNeedsFiles','MakeVideo=true 需要 WriteFiles=true。');
end
fontSize = 20; nTicks = 6; axisPad = 0.025;
if p.PlotResults
    [fig1,ax1] = paper_fixed_window('自由端 Z 坐标',p);
    hold(ax1,'on');
    plot(ax1,R.t,R.z0_vals(:,end),'LineWidth',2.5,'Color','#D95319');
    paper_style_axes(ax1,fontSize);
    xlabel(ax1,'$\tau$','Interpreter','latex','FontSize',fontSize);
    if strcmpi(p.PlotMode,'Relative')
        ylabel(ax1,'$\bar{z}^{(0)}(1,\tau)-c_z(\tau)$', ...
            'Interpreter','latex','FontSize',fontSize);
    else
        ylabel(ax1,'$\bar{z}^{(0)}(1,\tau)$', ...
            'Interpreter','latex','FontSize',fontSize);
    end
    [xl,xt] = adaptive_axis_ticks(R.t(:),nTicks,axisPad);
    [yl,yt] = adaptive_axis_ticks(R.z0_vals(:,end),nTicks,axisPad);
    set(ax1,'XLim',xl,'YLim',yl,'XTick',xt,'YTick',yt);
    title(ax1,'');
    pbaspect(ax1,[p.PlotBoxAspectRatio,1,1]);
    hold(ax1,'off');
    set(fig1,'Visible',p.FigureVisible);

    [shapeTimes,shapeRows] = paper_sample_rows(R.t,p.ShapeTimes);
    if ~isempty(shapeTimes)
        xp = R.x0_vals(shapeRows,:); zp = R.z0_vals(shapeRows,:);
        colors = turbo(numel(shapeTimes));
        [fig2,ax2] = paper_fixed_window('空间形变姿态叠加',p, ...
            p.ShapeFigureSize,p.ShapePlotBoxAspectRatio);
        hold(ax2,'on');
        paper_style_axes(ax2,fontSize);
        shapeHandles = gobjects(numel(shapeTimes),1);
        for k = 1:numel(shapeTimes)
            shapeHandles(k) = plot(ax2,xp(k,:),zp(k,:),'LineWidth',2.5,'Color',colors(k,:), ...
                'DisplayName',sprintf('$\\tau=%.15g$',shapeTimes(k)));
            plot(ax2,xp(k,1),zp(k,1),'o','Color',colors(k,:), ...
                'MarkerSize',7,'HandleVisibility','off');
            plot(ax2,xp(k,end),zp(k,end),'^','Color',colors(k,:), ...
                'MarkerSize',8,'MarkerFaceColor',colors(k,:),'HandleVisibility','off');
        end
        paper_shape_labels(ax2,p.PlotMode,fontSize);
        title(ax2,'');
        paper_fixed_shape_layout(ax2,shapeHandles,xp,zp,p,nTicks);
        set(fig2,'Visible',p.FigureVisible);
        if p.WriteFiles&&p.ExportShapeImage
            shapeImage=fullfile(p.OutputDir,sprintf('Rod_Shapes_square_%s.png',p.RunStamp));
            paper_assert_new_files({shapeImage});
            paper_export_square(fig2,ax2,shapeImage);
        end
    end
end

% 原程序即使不生成 MP4，也显示动画窗口的初始构型。
tauFrames = R.t(1):(p.PlaySpeed/p.FPS):R.t(end);
if isempty(tauFrames), tauFrames = R.t(1); end
if tauFrames(end) < R.t(end), tauFrames(end+1) = R.t(end); end
if numel(R.t) == 1
    xv = R.x0_vals; zv = R.z0_vals;
else
    xv = interp1(R.t,R.x0_vals,tauFrames,'pchip');
    zv = interp1(R.t,R.z0_vals,tauFrames,'pchip');
end
[figVideo,axv] = paper_fixed_window('杆件构型演化视频动画',p);
if ~p.PlotResults
    closeVideoFigure = onCleanup(@() paper_close_figure(figVideo)); 
end
hold(axv,'on');
paper_style_axes(axv,fontSize); paper_shape_labels(axv,p.PlotMode,fontSize);
[xl,yl,xt,yt] = paper_equal_limits(xv,zv,nTicks,axisPad,axisPad,p.PlotBoxAspectRatio);
axis(axv,'equal');
set(axv,'XLim',xl,'YLim',yl,'XTick',xt,'YTick',yt); title(axv,'');
pbaspect(axv,[p.PlotBoxAspectRatio,1,1]);
hRod = plot(axv,NaN,NaN,'LineWidth',2.5,'Color','#0072BD');
hRoot = plot(axv,NaN,NaN,'o','MarkerSize',8,'Color','#D95319','MarkerFaceColor','#D95319');
hTip = plot(axv,NaN,NaN,'^','MarkerSize',8,'Color','#EDB120','MarkerFaceColor','#EDB120');
hTime = text(axv,0.035,0.96,'','Units','normalized', ...
    'FontSize',fontSize,'FontWeight','normal','Interpreter','latex', ...
    'HorizontalAlignment','left','VerticalAlignment','top','Clipping','on');
set(figVideo,'Visible',p.FigureVisible);
if p.MakeVideo
    if ~exist(p.OutputDir,'dir'), mkdir(p.OutputDir); end
    videoName = fullfile(p.OutputDir,sprintf('Rod_Evolution_%s.mp4',p.RunStamp));
    paper_assert_new_files({videoName});
    videoObj = VideoWriter(videoName,'MPEG-4'); videoObj.FrameRate = p.FPS;
    open(videoObj);
    closeVideoFile = onCleanup(@() close(videoObj));
    fprintf('正在生成视频文件 %s ...\n',videoName);
    nFrames = numel(tauFrames);
else
    nFrames = 1;
end
for k = 1:nFrames
    set(hRod,'XData',xv(k,:),'YData',zv(k,:));
    set(hRoot,'XData',xv(k,1),'YData',zv(k,1));
    set(hTip,'XData',xv(k,end),'YData',zv(k,end));
    set(hTime,'String',sprintf('$\\tau=%.2f$',tauFrames(k)));
    drawnow;
    if p.MakeVideo, writeVideo(videoObj,getframe(figVideo)); end
end
if p.MakeVideo
    clear closeVideoFile;
    fprintf('视频生成完毕！保存至: %s\n',videoName);
end
end

function [times,rows] = paper_sample_rows(t,requestedTimes)
% 未到达的时间直接跳过；已到达但未包含在 tspan 中则明确报错。
t = t(:); requestedTimes = requestedTimes(:);
tol = 64*eps(max(1,max(abs(t))));
times = requestedTimes(requestedTimes >= t(1)-tol & requestedTimes <= t(end)+tol);
rows = zeros(numel(times),1);
for k = 1:numel(times)
    [distance,rows(k)] = min(abs(t-times(k)));
    if distance > tol
        error('beam:MissingOutputTime', ...
            '输出时刻 tau=%.15g 未包含在求解输出中，请将其加入输出 tspan。',times(k));
    end
end
end

function paper_assert_new_files(names)
for k = 1:numel(names)
    if exist(names{k},'file')
        error('beam:OutputExists', ...
            '为保护已有数据，不覆盖文件：%s。请更换 OutputDir 或稍后重新运行。',names{k});
    end
end
end

function paper_style_axes(ax,fontSize)
% 直接使用原生图框和刻度，不创建第二套右轴/上轴，也不手动画边框线。
box(ax,'on'); grid(ax,'off');
set(ax,'FontSize',fontSize,'TickLabelInterpreter','latex','LineWidth',1.2, ...
    'TickDir','in','XMinorTick','off','YMinorTick','off', ...
    'Color','w','XColor','k','YColor','k');
end

function [fig,ax]=paper_fixed_window(windowName,p,wh,aspect)
% 各类图窗使用预设像素尺寸和绘图区比例，不根据数据跨度改变窗口大小。
% 两参数调用用于末端位移图/视频；形态叠加图单独传入正方形窗口与绘图区。
if nargin<3,wh=p.FigureSize;end
if nargin<4,aspect=p.PlotBoxAspectRatio;end
fig=figure('Name',windowName,'NumberTitle','off','WindowStyle','normal','WindowState','normal', ...
    'Units','pixels','Position',[70,70,wh],'Resize','off','Visible','off', ...
    'Color','w','MenuBar','figure','ToolBar','none','InvertHardcopy','off');
% 固定留出纵轴标题、刻度和横轴标题所需空间；剩余空间按固定比例安排。
left=110; bottom=85; right=35; top=30;
available=wh-[left+right,bottom+top];
height=min(available(2),available(1)/aspect); width=aspect*height;
pos=[left+(available(1)-width)/2,bottom+(available(2)-height)/2,width,height];
ax=axes('Parent',fig,'Units','pixels','Position',pos, ...
    'PositionConstraint','innerposition','PlotBoxAspectRatio',[aspect,1,1]);
% 手动保存/打印时采用同一物理页面比例，不使用系统默认纸张尺寸。
screenDpi=get(groot,'ScreenPixelsPerInch'); paperSize=wh/screenDpi;
set(fig,'PaperUnits','inches','PaperSize',paperSize, ...
    'PaperPosition',[0,0,paperSize],'PaperPositionMode','manual');
end

function report=paper_fixed_shape_layout(ax,handles,x,z,p,nTicks)
% 固定图框、等跨度坐标和导出页面；不同时手动锁定 DAR 与 PBAR。
fig=ancestor(ax,'figure');
side=max(p.ShapeFigureSize);
set(fig,'WindowStyle','normal','WindowState','normal','Units','pixels', ...
    'Resize','off','Tag','BeamShapeSquare_v3');
fp=fig.Position;fig.Position=[fp(1:2),side,side];
set(ax,'Units','pixels','PositionConstraint','innerposition', ...
    'Position',[110,85,side-145,side-145]);
view(ax,2);
set(ax,'CameraViewAngleMode','auto','XScale','linear','YScale','linear');
[xl,yl,xt,yt]=paper_equal_limits(x,z,nTicks,[.025,.07],[.07,.07],1);
set(ax,'XLim',xl,'YLim',yl,'XTick',xt,'YTick',yt);
paper_square_limits(ax,nTicks);
for pass=1:4
    paper_fit_square_axes(fig,ax);
    paper_adaptive_legend(ax,handles,x,z,p.LegendFontSize,nTicks);
    % 必须在图例排版之后再次保证横纵跨度相等。
    paper_square_limits(ax,nTicks);
    drawnow;
    pos=ax.Position;inset=ax.TightInset;
    bounds=[pos(1:2)-inset(1:2),pos(1:2)+pos(3:4)+inset(3:4)];
    if all(bounds(1:2)>=4)&&all(bounds(3:4)<=side-4),break;end
end
% print 使用整张正方形页面，不对坐标轴单独裁边，也不改变宽高比。
set(fig,'PaperUnits','centimeters','PaperSize',[18,18], ...
    'PaperPosition',[0,0,18,18],'PaperPositionMode','manual','InvertHardcopy','off');
report=paper_square_check(fig,ax);
fprintf('[SquarePlot v3] 窗口宽/高=%.6f，图框宽/高=%.6f，横纵单位长度比=%.6f\n', ...
    report.canvas_ratio,report.frame_ratio,report.unit_ratio);
end

function paper_square_limits(ax,nTicks)
% X/Y 跨度相同且 PBAR=1，横纵单位长度自然相同；DAR 交给 MATLAB 计算。
xl=ax.XLim;yl=ax.YLim;
span=max(diff(xl),diff(yl));
xl=mean(xl)+[-.5,.5]*span;yl=mean(yl)+[-.5,.5]*span;
[~,xt]=adaptive_axis_ticks(xl,nTicks,0);
[~,yt]=adaptive_axis_ticks(yl,nTicks,0);
set(ax,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1,1,1], ...
    'PlotBoxAspectRatioMode','manual','XLim',xl,'YLim',yl, ...
    'XLimMode','manual','YLimMode','manual','XTick',xt,'YTick',yt);
end

function paper_fit_square_axes(fig,ax)
% 在固定画布内按实际文字占用留边，内图框始终为正方形。
drawnow;
wh=fig.Position(3:4);
margins=max([20,20,12,12],ceil(ax.TightInset+6));
available=wh-[sum(margins([1,3])),sum(margins([2,4]))];
side=floor(min(available));
if side<120,error('BeamGeneral:SquareCanvas','画布无法容纳坐标文字，请增大 ShapeFigureSize。');end
set(ax,'Position',[margins(1:2)+(available-side)/2,side,side]);
end

function report=paper_square_check(fig,ax)
% 校验实际布局；比例不满足时明确报错，不静默输出宽矩形。
drawnow;
fp=getpixelposition(fig);ap=getpixelposition(ax);
pb=ax.PlotBoxAspectRatio;da=ax.DataAspectRatio;
report=struct('canvas_ratio',fp(3)/fp(4),'frame_ratio',ap(3)/ap(4), ...
    'unit_ratio',(ap(3)/diff(ax.XLim))/(ap(4)/diff(ax.YLim)));
ratios=[report.canvas_ratio,report.frame_ratio,report.unit_ratio, ...
    pb(1)/pb(2),da(1)/da(2)];
if any(~isfinite(ratios))||any(abs(ratios-1)>1e-6)
    error('BeamGeneral:SquareCheck', ...
        '正方形校验失败：窗口=%.6g，图框=%.6g，单位长度比=%.6g。', ...
        report.canvas_ratio,report.frame_ratio,report.unit_ratio);
end
end

function paper_export_square(fig,ax,filename)
paper_square_check(fig,ax);
[folder,~,ext]=fileparts(filename);
if isempty(ext),filename=[filename,'.png'];
elseif ~strcmpi(ext,'.png'),error('BeamGeneral:ImageFormat','请使用 .png 文件名。');end
if ~isempty(folder)&&~exist(folder,'dir'),mkdir(folder);end
set(fig,'PaperUnits','centimeters','PaperSize',[18,18], ...
    'PaperPosition',[0,0,18,18],'PaperPositionMode','manual');
print(fig,filename,'-dpng','-r600');
info=imfinfo(filename);
if abs(info.Width-info.Height)>1
    error('BeamGeneral:SquareExport','导出图像不是正方形：%d x %d。',info.Width,info.Height);
end
paper_square_check(fig,ax);
fprintf('[SquarePlot v3] 已导出 %d x %d 像素：%s\n',info.Width,info.Height,filename);
end


function paper_adaptive_legend(ax,handles,x,z,fontSize,nTicks)
% 先尝试不同列数及框内位置；仅在没有空位时扩展数据范围，图窗尺寸不变。
% 用线段/矩形相交检测，避免仅检查离散节点而漏掉跨过图例的曲线。
labels=cell(numel(handles),1);
for j=1:numel(handles),labels{j}=handles(j).DisplayName;end
lgd=legend(ax,handles,labels,'Interpreter','latex','FontSize',fontSize, ...
    'Box','on','Color','w','TextColor','k','LineWidth',.8, ...
    'AutoUpdate','off','Orientation','vertical','Units','pixels','Location','northeast');
drawnow;
box=ax.Position; boxSize=box(3:4); xl=ax.XLim; yl=ax.YLim;
pixelsPerPoint=get(groot,'ScreenPixelsPerInch')/72;
tickPixels=ax.TickLength(1)*max(boxSize);
inset=ceil(tickPixels+max(6,.25*ax.FontSize*pixelsPerPoint)+ ...
    .5*(ax.LineWidth+lgd.LineWidth)*pixelsPerPoint);
clearance=ceil(max(10,.45*ax.FontSize*pixelsPerPoint));
inset=max(inset,clearance+2);
segments=paper_pixel_segments(x,z,xl,yl,boxSize);
bestScore=inf; bestArea=inf; best=struct();
layouts=struct('columns',{},'size',{});

% 优先保留用户设定的字号；只有所有列数都放不下时才逐级减小。
for currentFont=fontSize*[1,.9,.8,.7]
    set(lgd,'FontSize',currentFont);
    for columns=1:min(numel(handles),8)
        set(lgd,'NumColumns',columns,'Location','northeast');drawnow;
        legendPos=lgd.Position; legendSize=ceil(legendPos(3:4));
        if any(legendSize>boxSize-2*inset),continue;end
        layouts(end+1)=struct('columns',columns,'size',legendSize); %#ok<AGROW>
        last=boxSize-inset-legendSize;
        [gx,gy]=meshgrid(linspace(inset,last(1),5),linspace(inset,last(2),5));
        positions=unique([inset,inset;last;inset,last(2);last(1),inset; ...
            gx(:),gy(:)],'rows','stable');
        for k=1:size(positions,1)
            rect=[positions(k,:),legendSize];
            protected=[rect(1:2)-clearance,rect(3:4)+2*clearance];
            score=paper_rectangle_overlap(segments,protected);
            area=prod(legendSize);
            if score<bestScore || (score==bestScore && area<bestArea)
                bestScore=score;bestArea=area;
                best=struct('columns',columns,'rect',rect,'font',currentFont);
            end
        end
    end
    if ~isempty(layouts),break;end
end
if isempty(layouts)
    error('BeamGeneral:LegendSize', ...
        '图例条目过多或过长，无法放入当前图框。请减少 ShapeTimes 或增大 ShapeFigureSize。');
end

if bestScore>0
    % 比较四边留出空白带所需的范围增量，选择对数据缩放影响最小的方案。
    % X/Y 跨度同步按同一倍数扩展，仍满足等单位长度和固定图框比例。
    finite=isfinite(x)&isfinite(z);
    bounds=[min(x(finite)),max(x(finite)),min(z(finite)),max(z(finite))];
    xSpan=diff(xl);zSpan=diff(yl);bestScale=inf;bestArea=inf;
    for j=1:numel(layouts)
        sz=layouts(j).size;
        fraction=(inset+sz+clearance)./boxSize;
        if any(fraction>=1),continue;end
        scales=[(xl(2)-bounds(1))/(xSpan*(1-fraction(1))), ...
            (bounds(2)-xl(1))/(xSpan*(1-fraction(1))), ...
            (yl(2)-bounds(3))/(zSpan*(1-fraction(2))), ...
            (bounds(4)-yl(1))/(zSpan*(1-fraction(2)))];
        scales=max(1,scales)*(1+1e-9);
        for side=1:4
            scale=scales(side);area=prod(sz);
            if scale<bestScale || (scale==bestScale && area<bestArea)
                bestScale=scale;bestArea=area;
                best.columns=layouts(j).columns;
                best.size=sz;best.side=side;
            end
        end
    end
    newX=mean(xl)+[-.5,.5]*xSpan*bestScale;
    newZ=mean(yl)+[-.5,.5]*zSpan*bestScale;
    sz=best.size;
    switch best.side
        case 1 % 左侧空白带
            newX=[xl(2)-xSpan*bestScale,xl(2)];xy=[inset,inset];
        case 2 % 右侧空白带
            newX=[xl(1),xl(1)+xSpan*bestScale];xy=[boxSize(1)-inset-sz(1),inset];
        case 3 % 下侧空白带
            newZ=[yl(2)-zSpan*bestScale,yl(2)];xy=[inset,inset];
        case 4 % 上侧空白带
            newZ=[yl(1),yl(1)+zSpan*bestScale];xy=[inset,boxSize(2)-inset-sz(2)];
    end
    [~,xt]=adaptive_axis_ticks(newX,nTicks,0);
    [~,yt]=adaptive_axis_ticks(newZ,nTicks,0);
    set(ax,'XLim',newX,'YLim',newZ,'XTick',xt,'YTick',yt);
    best.rect=[xy,sz];
end
set(lgd,'FontSize',best.font,'NumColumns',best.columns,'Location','none');
drawnow;
lgd.Position=[box(1:2)+best.rect(1:2),best.rect(3:4)];
end

function segments=paper_pixel_segments(x,z,xl,yl,boxSize)
% 每行是一条构型；分别提取相邻节点，禁止连接不同时刻的两条曲线。
px=(x-xl(1))*boxSize(1)/diff(xl);
pz=(z-yl(1))*boxSize(2)/diff(yl);
x1=px(:,1:end-1);x2=px(:,2:end);z1=pz(:,1:end-1);z2=pz(:,2:end);
segments=[x1(:),z1(:),x2(:),z2(:)];
segments=segments(all(isfinite(segments),2),:);
end

function score=paper_rectangle_overlap(segments,rect)
% 参数化线段裁切；同时计入相交段数和覆盖长度，零分表示完全没有相交。
count=size(segments,1);enter=zeros(count,1);leave=ones(count,1);valid=true(count,1);
for direction=1:2
    first=segments(:,direction);delta=segments(:,direction+2)-first;
    lower=rect(direction);upper=lower+rect(direction+2);moving=delta~=0;
    valid=valid & (moving | (first>=lower & first<=upper));
    a=-inf(count,1);b=inf(count,1);
    t1=(lower-first(moving))./delta(moving);t2=(upper-first(moving))./delta(moving);
    a(moving)=min(t1,t2);b(moving)=max(t1,t2);
    enter=max(enter,a);leave=min(leave,b);
end
hit=valid & enter<=leave;
lengths=hypot(segments(:,3)-segments(:,1),segments(:,4)-segments(:,2));
score=nnz(hit)+sum(lengths(hit).*max(0,leave(hit)-enter(hit)));
end

function paper_shape_labels(ax,plotMode,fontSize)
if strcmpi(plotMode,'Relative')
    xlabel(ax,'$\bar{x}^{(0)}-c_x$','Interpreter','latex','FontSize',fontSize);
    ylabel(ax,'$\bar{z}^{(0)}-c_z$','Interpreter','latex','FontSize',fontSize);
else
    xlabel(ax,'$\bar{x}^{(0)}$','Interpreter','latex','FontSize',fontSize);
    ylabel(ax,'$\bar{z}^{(0)}$','Interpreter','latex','FontSize',fontSize);
end
end

function paper_close_figure(fig)
if isgraphics(fig), close(fig); end
end

function [xl,yl,xt,yt]=paper_equal_limits(x,z,nTicks,xPad,zPad,boxAspect)
% 固定绘图区宽高比，并让 X/Z 每单位对应相同像素长度。
% 图框横向较宽时，X 的数据跨度按相同比例增大，不拉伸梁的形状。
[xl,~]=adaptive_axis_ticks(x(:),nTicks,xPad);
[yl,~]=adaptive_axis_ticks(z(:),nTicks,zPad);
zSpan=max(diff(yl),diff(xl)/boxAspect);
xl=mean(xl)+[-.5,.5]*(boxAspect*zSpan);yl=mean(yl)+[-.5,.5]*zSpan;
[~,xt]=adaptive_axis_ticks(xl,nTicks,0);[~,yt]=adaptive_axis_ticks(yl,nTicks,0);
end

function [axisLim,tickValues] = adaptive_axis_ticks(data,targetTicks,padFraction)
% 原程序的自适应坐标范围与 5--6 个主刻度算法。
data = data(isfinite(data));
if isempty(data)
    axisLim = [-1,1]; tickValues = linspace(-1,1,targetTicks); return;
end
dmin = min(data); dmax = max(data); dataSpan = dmax-dmin;
dataScale = max([abs(dmin),abs(dmax),1]);
if dataSpan <= 1e-12*dataScale
    halfWidth = max(0.025*dataScale,1e-4);
    dmin = dmin-halfWidth; dmax = dmax+halfWidth; dataSpan = dmax-dmin;
end
if isscalar(padFraction)
    padLower = padFraction; padUpper = padFraction;
elseif numel(padFraction) == 2
    padLower = padFraction(1); padUpper = padFraction(2);
else
    error('padFraction must be a scalar or a two-element vector [lowerPad, upperPad].');
end
axisLim = [dmin-padLower*dataSpan,dmax+padUpper*dataSpan];
rawStep = dataSpan/max(targetTicks-1,1); exponent = floor(log10(rawStep));
base = 10^exponent; candidates = base*[0.2,0.25,0.5,1,2,2.5,5,10];
bestScore = inf; bestTicks = [];
for step = candidates
    firstTick = ceil((axisLim(1)-1e-12*abs(axisLim(1)))/step)*step;
    lastTick = floor((axisLim(2)+1e-12*abs(axisLim(2)))/step)*step;
    if lastTick < firstTick, continue; end
    ticks = firstTick:step:lastTick; nTick = numel(ticks);
    if nTick < 3 || nTick > 8, continue; end
    if nTick == 5 || nTick == 6
        countPenalty = 0;
    else
        countPenalty = abs(nTick-5.5);
    end
    stepPenalty = 0.05*abs(log(step/rawStep)); score = countPenalty+stepPenalty;
    if score < bestScore, bestScore = score; bestTicks = ticks; end
end
if isempty(bestTicks), bestTicks = linspace(axisLim(1),axisLim(2),targetTicks); end
tickValues = bestTicks; tickValues(abs(tickValues)<1e-12) = 0;
end
