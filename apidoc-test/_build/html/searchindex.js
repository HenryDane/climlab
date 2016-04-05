Search.setIndex({envversion:46,filenames:["climlab","climlab.convection","climlab.domain","climlab.dynamics","climlab.model","climlab.process","climlab.radiation","climlab.solar","climlab.surface","climlab.tests","climlab.utils","index","modules","setup"],objects:{"":{climlab:[0,0,0,"-"]},"climlab.convection":{convadj:[1,0,0,"-"]},"climlab.convection.convadj":{Akamaev_adjustment:[1,2,1,""],Akamaev_adjustment_multidim:[1,2,1,""],ConvectiveAdjustment:[1,3,1,""],convective_adjustment_direct:[1,2,1,""]},"climlab.convection.convadj.ConvectiveAdjustment":{adj_lapse_rate:[1,4,1,""],compute:[1,1,1,""]},"climlab.domain":{axis:[2,0,0,"-"],domain:[2,0,0,"-"],field:[2,0,0,"-"]},"climlab.domain.axis":{Axis:[2,3,1,""]},"climlab.domain.domain":{Atmosphere:[2,3,1,""],Ocean:[2,3,1,""],SlabAtmosphere:[2,3,1,""],SlabOcean:[2,3,1,""],box_model_domain:[2,2,1,""],make_slabatm_axis:[2,2,1,""],make_slabocean_axis:[2,2,1,""],single_column:[2,2,1,""],zonal_mean_column:[2,2,1,""],zonal_mean_surface:[2,2,1,""]},"climlab.domain.domain.Atmosphere":{set_heat_capacity:[2,1,1,""]},"climlab.domain.domain.Ocean":{set_heat_capacity:[2,1,1,""]},"climlab.domain.field":{Field:[2,3,1,""],global_mean:[2,2,1,""]},"climlab.dynamics":{budyko_transport:[3,0,0,"-"],diffusion:[3,0,0,"-"]},"climlab.dynamics.budyko_transport":{BudykoTransport:[3,3,1,""]},"climlab.dynamics.budyko_transport.BudykoTransport":{b:[3,4,1,""]},"climlab.dynamics.diffusion":{Diffusion:[3,3,1,""],MeridionalDiffusion:[3,3,1,""]},"climlab.model":{column:[4,0,0,"-"],ebm:[4,0,0,"-"]},"climlab.model.column":{BandRCModel:[4,3,1,""],GreyRadiationModel:[4,3,1,""],RadiativeConvectiveModel:[4,3,1,""],compute_layer_absorptivity:[4,2,1,""],initial_state:[4,2,1,""]},"climlab.model.column.GreyRadiationModel":{do_diagnostics:[4,1,1,""],initial_state:[4,1,1,""]},"climlab.model.ebm":{EBM:[4,3,1,""],EBM_annual:[4,3,1,""],EBM_seasonal:[4,3,1,""]},"climlab.model.ebm.EBM":{diffusive_heat_transport:[4,1,1,""],global_mean_temperature:[4,1,1,""],heat_transport:[4,1,1,""],heat_transport_convergence:[4,1,1,""],inferred_heat_transport:[4,1,1,""]},"climlab.process":{diagnostic:[5,0,0,"-"],energy_budget:[5,0,0,"-"],implicit:[5,0,0,"-"],process:[5,0,0,"-"],time_dependent_process:[5,0,0,"-"]},"climlab.process.diagnostic":{DiagnosticProcess:[5,3,1,""]},"climlab.process.diagnostic.DiagnosticProcess":{compute:[5,1,1,""]},"climlab.process.energy_budget":{EnergyBudget:[5,3,1,""],ExternalEnergySource:[5,3,1,""]},"climlab.process.energy_budget.EnergyBudget":{compute:[5,1,1,""]},"climlab.process.implicit":{ImplicitProcess:[5,3,1,""]},"climlab.process.implicit.ImplicitProcess":{compute:[5,1,1,""]},"climlab.process.process":{Process:[5,3,1,""],get_axes:[5,2,1,""],process_like:[5,2,1,""]},"climlab.process.process.Process":{add_subprocess:[5,1,1,""],add_subprocesses:[5,1,1,""],depth:[5,4,1,""],depth_bounds:[5,4,1,""],lat:[5,4,1,""],lat_bounds:[5,4,1,""],lev:[5,4,1,""],lev_bounds:[5,4,1,""],lon:[5,4,1,""],lon_bounds:[5,4,1,""],remove_subprocess:[5,1,1,""],set_state:[5,1,1,""]},"climlab.process.time_dependent_process":{TimeDependentProcess:[5,3,1,""]},"climlab.process.time_dependent_process.TimeDependentProcess":{compute:[5,1,1,""],compute_diagnostics:[5,1,1,""],integrate_converge:[5,1,1,""],integrate_days:[5,1,1,""],integrate_years:[5,1,1,""],set_timestep:[5,1,1,""],step_forward:[5,1,1,""]},"climlab.radiation":{AplusBT:[6,0,0,"-"],cloud:[6,0,0,"-"],insolation:[6,0,0,"-"],nband:[6,0,0,"-"],radiation:[6,0,0,"-"],transmissivity:[6,0,0,"-"],water_vapor:[6,0,0,"-"]},"climlab.radiation.AplusBT":{AplusBT:[6,3,1,""]},"climlab.radiation.AplusBT.AplusBT":{A:[6,4,1,""],B:[6,4,1,""]},"climlab.radiation.cloud":{Reflection:[6,2,1,""],compute_beta:[6,2,1,""],compute_eps:[6,2,1,""],compute_tauN:[6,2,1,""]},"climlab.radiation.insolation":{AnnualMeanInsolation:[6,3,1,""],DailyInsolation:[6,3,1,""],FixedInsolation:[6,3,1,""],P2Insolation:[6,3,1,""]},"climlab.radiation.insolation.AnnualMeanInsolation":{orb:[6,4,1,""]},"climlab.radiation.insolation.P2Insolation":{s2:[6,4,1,""]},"climlab.radiation.nband":{FourBandLW:[6,3,1,""],FourBandSW:[6,3,1,""],NbandRadiation:[6,3,1,""],SPEEDY_band_fraction:[6,2,1,""],ThreeBandSW:[6,3,1,""]},"climlab.radiation.nband.FourBandSW":{emissivity:[6,4,1,""]},"climlab.radiation.nband.NbandRadiation":{band_fraction:[6,4,1,""]},"climlab.radiation.nband.ThreeBandSW":{emissivity:[6,4,1,""]},"climlab.radiation.radiation":{Radiation:[6,3,1,""],RadiationLW:[6,3,1,""],RadiationSW:[6,3,1,""]},"climlab.radiation.transmissivity":{Transmissivity:[6,3,1,""],compute_T_vectorized:[6,2,1,""],tril:[6,2,1,""]},"climlab.radiation.transmissivity.Transmissivity":{flux_down:[6,1,1,""],flux_reflected_up:[6,1,1,""],flux_up:[6,1,1,""]},"climlab.radiation.water_vapor":{FixedRelativeHumidity:[6,3,1,""],ManabeWaterVapor:[6,3,1,""]},"climlab.solar":{insolation:[7,0,0,"-"],orbital:[7,0,0,"-"],orbital_cycles:[7,0,0,"-"]},"climlab.solar.insolation":{daily_insolation:[7,2,1,""],solar_longitude:[7,2,1,""]},"climlab.solar.orbital":{LongOrbitalTable:[7,3,1,""],OrbitalTable:[7,3,1,""]},"climlab.solar.orbital.OrbitalTable":{lookup_parameters:[7,1,1,""]},"climlab.solar.orbital_cycles":{OrbitalCycles:[7,3,1,""]},"climlab.surface":{albedo:[8,0,0,"-"],surface_radiation:[8,0,0,"-"],turbulent:[8,0,0,"-"]},"climlab.surface.albedo":{ConstantAlbedo:[8,3,1,""],Iceline:[8,3,1,""],P2Albedo:[8,3,1,""],StepFunctionAlbedo:[8,3,1,""]},"climlab.surface.albedo.ConstantAlbedo":{albedo:[8,4,1,""]},"climlab.surface.albedo.Iceline":{compute:[8,1,1,""],find_icelines:[8,1,1,""]},"climlab.surface.albedo.P2Albedo":{a0:[8,4,1,""],a2:[8,4,1,""]},"climlab.surface.albedo.StepFunctionAlbedo":{compute:[8,1,1,""]},"climlab.surface.surface_radiation":{SurfaceRadiation:[8,3,1,""]},"climlab.surface.surface_radiation.SurfaceRadiation":{emission:[8,1,1,""],reflected_flux:[8,1,1,""]},"climlab.surface.turbulent":{LatentHeatFlux:[8,3,1,""],SensibleHeatFlux:[8,3,1,""],SurfaceFlux:[8,3,1,""]},"climlab.surface.turbulent.LatentHeatFlux":{compute_flux:[8,1,1,""]},"climlab.surface.turbulent.SensibleHeatFlux":{compute_flux:[8,1,1,""]},"climlab.utils":{constants:[10,0,0,"-"],heat_capacity:[10,0,0,"-"],legendre:[10,0,0,"-"],thermo:[10,0,0,"-"],walk:[10,0,0,"-"]},"climlab.utils.heat_capacity":{atmosphere:[10,2,1,""],ocean:[10,2,1,""],slab_ocean:[10,2,1,""]},"climlab.utils.legendre":{P0:[10,2,1,""],P10:[10,2,1,""],P10prime:[10,2,1,""],P12:[10,2,1,""],P12prime:[10,2,1,""],P14:[10,2,1,""],P14prime:[10,2,1,""],P16:[10,2,1,""],P18:[10,2,1,""],P1:[10,2,1,""],P1prime:[10,2,1,""],P20:[10,2,1,""],P22:[10,2,1,""],P24:[10,2,1,""],P26:[10,2,1,""],P28:[10,2,1,""],P2:[10,2,1,""],P2prime:[10,2,1,""],P3:[10,2,1,""],P3prime:[10,2,1,""],P4:[10,2,1,""],P4prime:[10,2,1,""],P5:[10,2,1,""],P6:[10,2,1,""],P6prime:[10,2,1,""],P8:[10,2,1,""],P8prime:[10,2,1,""],Pn:[10,2,1,""],Pnprime:[10,2,1,""]},"climlab.utils.thermo":{EIS:[10,2,1,""],Planck_frequency:[10,2,1,""],Planck_wavelength:[10,2,1,""],Planck_wavenumber:[10,2,1,""],T:[10,2,1,""],blackbody_emission:[10,2,1,""],clausius_clapeyron:[10,2,1,""],estimated_inversion_strength:[10,2,1,""],potential_temperature:[10,2,1,""],pseudoadiabat:[10,2,1,""],qsat:[10,2,1,""],temperature_from_potential:[10,2,1,""],theta:[10,2,1,""]},"climlab.utils.walk":{process_tree:[10,2,1,""],walk_processes:[10,2,1,""]},climlab:{convection:[1,0,0,"-"],domain:[2,0,0,"-"],dynamics:[3,0,0,"-"],model:[4,0,0,"-"],process:[5,0,0,"-"],radiation:[6,0,0,"-"],solar:[7,0,0,"-"],surface:[8,0,0,"-"],tests:[9,0,0,"-"],utils:[10,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","method","Python method"],"2":["py","function","Python function"],"3":["py","class","Python class"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:module","1":"py:method","2":"py:function","3":"py:class","4":"py:attribute"},terms:{"30degc":10,"45n":4,"65n":7,"_domain":2,"_insol":6,"abstract":2,"class":[1,2,3,4,5,6,7,8,10],"default":[1,2,3,4,5,6,7],"final":6,"function":10,"import":[3,4,6,7],"long":[4,7],"new":[2,4,5],"return":[1,2,4,5,6,7,10],"short":10,"super":10,"true":[5,7,10],"try":7,abov:[6,7],abs_coeff:4,absolut:10,absorb:6,absorbed_tot:6,absorber_vmr:6,absorpt:[4,6],absorption_cross_sect:6,access:7,accord:7,account:5,accur:10,actual:6,add:[4,5,7],add_subprocess:[4,5,7],adj_lapse_r:[1,4],adjust:1,after:7,aim_v23:6,air:[6,10],akamaev_adjust:1,akamaev_adjustment_multidim:1,albani:[4,10],albedo:[],albedo_sfc:[4,6,8],all:[4,5,6,7,10],alon:3,along:[6,10],also:[2,5,6],alwai:7,angl:7,angular:7,ani:7,annualmeaninsol:6,aplusbt:[],aplusbt_co2:[],append:10,approxim:7,area:10,argument:[1,5,7],arrai:[2,3,6,7],ascent:10,asd:7,associ:2,assum:[6,10],atm:[2,6,10],atmopsher:1,atmospher:[2,6,7,10],attempt:7,attribut:[2,6],august:7,author:[6,7],automat:7,avail:7,averag:[4,7],axi:[],axial:7,axis_typ:[2,4],balanc:4,band:[3,6],band_fract:6,bandrcmodel:4,base:[1,2,3,4,5,6,7,8,10],basic:6,beam:6,becaus:[6,7],befor:7,behav:2,berger:7,beta:1,between:[3,6,10],black:[],blackbodi:10,blackbody_emiss:10,bodi:[],bolton:10,boltzmann:[],both:[7,10],bottom:6,bound:2,boundari:6,box:[2,10],box_model_domain:2,bretherton:10,brian:[4,6,10],broadcast:6,brose:[4,10],budget:5,budyko:3,budyko_transport:[],budykotransport:3,bunch:6,calcul:[2,3,4,6,7,10],caldeira:[],calendar:[5,7],call:[2,5],can:[2,3,5,6,7],capac:[1,10],chang:[5,7,10],characterist:6,claim:10,clausius_clapeyron:10,climat:[7,10],climtrad:[],cloud:[],co2:[],code:[4,6,7,10],coeffici:4,col:4,collect:10,column:[],common:10,complet:10,comput:[1,4,5,6,7,8,10],compute_absorpt:[],compute_beta:6,compute_diagnost:5,compute_emiss:[],compute_ep:6,compute_flux:8,compute_layer_absorpt:4,compute_optical_path:[],compute_t_vector:6,compute_taun:6,concentr:[],conserv:7,consid:[],constant:[],constantalbedo:8,construct:6,contain:[5,7],contribut:[],convadj:[],convect:[],convective_adjustment_direct:1,convectiveadjust:1,conveni:[2,4,10],convent:7,converg:[4,5,6],convert:10,copi:[3,6],correspond:5,coszen:6,coupl:5,cours:10,creat:[2,4,6,10],creation:5,crit:5,criteria:5,cross:6,cumul:6,current:[5,7],custom:2,cycl:7,dai:[5,7],daili:7,daily_insol:7,dailyinsol:[4,6],data:7,date:5,day_typ:7,days_per_year:7,defin:7,definit:10,deg2rad:3,degc:[3,10],degre:[1,4,7],densiti:10,depend:[2,5],depth:[2,5,10],depth_bound:5,deriv:10,describ:10,design:7,develop:[4,10],diagnost:[],diagnosticprocess:[5,6,8],diagon:6,dictionari:[1,5,6,7,10],diff:6,differ:[3,5],differenti:7,diffus:[],diffusion_axi:3,diffusive_heat_transport:4,dimens:[6,7],dimension:[3,6],dimensionless:[7,10],discret:[5,6],distanc:7,divid:6,do_diagnost:4,doe:6,domain:[],don:5,down:6,downward:1,downwel:6,due:[],dynam:[],each:[5,6,7],earth:7,ebm:[],ebm_annu:4,ebm_season:[4,7],ecc:[6,7],eccentr:7,edown:6,edu:[4,7,10],eisenman:7,either:5,element:6,emiss:[6,8],emissivity_sfc:6,emit_sfc:6,enddo:6,energi:[4,5],energy_budget:[],energybudget:[3,4,5,6,8],env:10,eps3:6,epslw:6,equal:7,equat:[3,5,6,7],equinox:7,equip:7,estim:[7,10],estimated_inversion_strength:10,eup:6,evenli:2,everi:[2,5,6],exactli:[2,7],exampl:[3,4,6,7],exit:5,experi:7,explicit:5,express:1,externalenergysourc:5,fals:[3,10],faster:3,fband:6,feedback:7,field:[],file:7,find_icelin:8,first:[7,10],fix:[4,5],fixedinsol:6,fixedrelativehumid:6,flag:3,flux:[5,6,10],flux_components_bottom:[],flux_components_top:[],flux_down:6,flux_from_sfc:6,flux_from_spac:6,flux_reflected_up:6,flux_to_sfc:6,flux_to_spac:6,flux_up:6,fluxdown:6,fluxdowntop:6,fluxupbottom:6,follow:[6,7,10],formula:[6,10],fortran:6,forward:5,found:7,four:6,fourbandlw:6,fourbandsw:6,fraction:6,frequenc:10,from:[2,3,4,5,6,7,10],fromspac:6,fsw:7,ftp:7,fulli:7,further:[],gener:[5,7,10],get:6,get_ax:5,give:[1,6],given:[3,5,6,7,10],global:[2,3,4],global_mean:2,global_mean_temperatur:4,gov:7,grei:4,greyga:6,greyradiationmodel:[3,4],grid:[2,4,5,6,10],handl:[7,10],harvard:7,have:7,heat:[1,3,4,5,10],heat_capac:[],heat_transport:4,heat_transport_converg:4,heating_r:5,here:[3,4,6],hold:7,hpa:[1,10],http:7,humid:10,huyber:7,ian:7,icelin:8,ident:[5,7,10],ignoreflag:10,imbal:4,imcc:7,implement:[3,4,6],implicit:[],implicitprocess:[3,5],includ:[1,7],incom:7,increas:1,index:11,infer:4,inferred_heat_transport:4,initi:3,initial_st:4,input:[1,3,5,6,7,10],insol:[],insola:7,instantan:[4,10],instantli:1,integr:[4,5,7],integrate_converg:5,integrate_dai:5,integrate_year:[3,5],interfac:6,interpol:7,interv:10,invers:10,invok:[6,7],isn:7,item:5,iter:5,jan:7,januari:7,journal:7,jtemp:6,just:[3,5],kast:[],kei:5,kelvin:[6,10],kepler:7,kreuzer:3,kwarg:[1,2,3,4,5,6,8],kyear:7,kyear_start:7,kyear_stop:7,la2004:7,laboratori:10,lambda:7,laps:1,lapser:1,laskar:7,last:[6,7],lat:[2,3,4,5,7],lat_bound:5,latentheatflux:8,latitud:[2,4,7],law:[7,10],layer:[4,6],legendr:[],length:[3,7],lev:[2,4,5],lev_bound:5,level:[2,6],like:[2,5],linalg:3,linear:[6,7],linearli:7,linspac:7,list:[2,10],load:7,local:[3,7,10],lon:[2,5],lon_bound:5,long_peri:[6,7],longitud:7,longorbitalt:7,longwav:6,look:7,lookup_paramet:7,loop:5,lot:10,loutr:7,lower:6,lwtemp1:6,lwtemp2:6,make:[3,5,6],make_slabatm_axi:2,make_slabocean_axi:2,manabewatervapor:6,mar:6,march:7,mass:6,matlab:7,matplotlib:3,matrix:6,mean:[1,2,3,4],measur:7,member:7,meridion:3,meridionaldiffus:3,meter:10,method:[2,4,7,10],million:7,minimum:6,mitgcm:6,mix:6,model:[],modifi:[5,7],momentum:7,monthli:10,moritz:3,multi:6,multipl:[6,7],myear:7,n_k:1,name:[5,10],nband:[],nbandradi:6,ncdc:7,ndarrai:2,necessari:[],need:6,neutral:1,noaa:7,none:[1,2,3,4,5,6,7,8],north:4,note:[4,7],now:[4,6],num_it:5,num_lat:[2,3,4,5],num_lev:[2,4],num_level:5,num_point:2,num_steps_per_year:5,number:[2,5,6],numpi:[2,3,6],object:[2,4,5,6],obliqu:[6,7],obtain:7,occur:[3,7],ocean:[2,10],olr:6,onli:[1,3,7],onlin:7,oop:5,oper:3,option:7,orb:[6,7],orbit91:7,orbit:[],orbital_cycl:[],orbital_year_factor:7,orbitalcycl:7,orbitalt:7,orient:4,origin:7,other:6,otherwis:1,out:10,outgo:[],output:[7,10],over:7,overli:2,overrid:[],p10:10,p10prime:10,p12:10,p12prime:10,p14:10,p14prime:10,p16:10,p18:10,p1prime:10,p20:10,p22:10,p24:10,p26:10,p28:10,p2albedo:8,p2insol:6,p2prime:10,p3prime:10,p4prime:10,p6prime:10,p8prime:10,page:11,paleo:7,paper:10,param:[3,4,5,7],paramet:[3,7],parcel:10,parent:[3,5,6,7],part:[6,10],partit:6,pass:[2,5,6],past:7,per:[1,5,6],pergammon:10,perihelion:7,peter:7,phy_radiat:6,physic:10,pierrehumbert:10,pkg:6,planck:10,planck_frequ:10,planck_wavelength:10,planck_wavenumb:10,planetari:10,plot:3,plt:3,plu:6,pnprime:10,point:[2,3,4],pole:4,polynomi:10,polyomi:10,port:7,posit:[1,6],potenti:10,potential_temperatur:10,ppm:[],precess:7,prescrib:4,present:7,press:10,pressur:[1,2,10],principl:10,print:[2,3,6,7],proc:[5,10],procdict:5,process:[],process_lik:5,process_or_domain:5,process_tre:10,product:6,provid:[6,7],pseudoadiabat:10,pub:7,purpos:7,pyplot:3,python:6,qsat:10,qstrat:6,quantiti:[2,5],quaternari:7,radi:[4,6],radiat:[],radiationlw:6,radiationsw:6,radiative_h:[],radiativeconvectivemodel:4,radit:10,rain:10,rate:[1,5,10],rather:3,ratio:6,raymond:10,read:7,readm:7,recurs:10,refer:7,referenc:7,reflect:6,reflected_flux:8,rel:7,relat:7,relative_humid:6,remain:3,remot:7,remov:5,remove_subprocess:5,replac:4,repres:[4,7],represent:10,reproduc:6,result:6,review:[7,10],roger:10,rose:[4,10],routin:[7,10],run:7,s65:7,s_k:1,same:[6,7],satur:10,scalar:[6,7],schwarschild:6,scienc:7,scipi:3,search:11,season:4,second:[5,7,10],section:[6,7],see:[],segment_length_year:7,self:6,semi:[],sensibleheatflux:8,set:[3,4,5,7],set_heat_capac:2,set_stat:5,set_timestep:5,setup:[],sever:10,sfc:[2,3,6],shortwav:4,should:[3,5,6,7],show:3,simpl:[2,6],simplest:6,sin:3,singl:[1,2,4,5,6],single_column:[2,6],sink:5,situ:10,size:7,slab:[2,10],slab_ocean:10,slabatmospher:2,slabocean:2,slope:10,solar:[],solar_longitud:7,solstic:7,solut:[5,7],solv:[3,5,6],solve_band:3,solver:3,some:10,sourc:[1,2,3,4,5,6,7,8,10],space:[2,6,10],spatial:[3,5],specif:[7,10],specifi:[1,5,6,7],spectral:6,spectrum:6,speedi:6,speedy_band_fract:6,split_channel:[],spring:7,stagger:4,stand:3,start:[6,7],state:[1,2,3,5,6,7],stefan:10,step:5,step_forward:[3,5,6],stepfunctionalbedo:[7,8],stommelbox:[],store:6,str:5,stream:6,strength:10,strictli:5,string:10,sub:10,subclass:7,submodul:[],subproces:5,subprocess:[3,4,5],subset:3,sum:6,summer:7,sun:7,support:10,sure:5,surfac:[],surface_radi:[],surfaceflux:8,surfaceradi:8,t700:10,t_init_0:[],t_init_p2:[],t_k:1,take:6,tatm:[3,6],tau0:6,tau1:6,tau2:6,tau:6,taun:6,tdown:6,temp:10,temperatur:[1,3,4,6,10],temperature_from_potenti:10,tendenc:[5,6],term:7,test:[],than:3,them:10,thermal:3,thermo:[],thermodynam:10,theta:[1,10],theta_k:1,thi:[2,4,5,7],thousand:7,three:7,threebandsw:6,through:[5,7],thu:[4,6],tilt:7,time:[3,4,5,7],time_dependent_process:[],time_typ:5,timedependentprocess:[1,4,5],timeseri:7,timestep:[4,5],toa:[2,4],top:[6,7,10],top_proc:10,topdown:[5,10],topnam:10,total:6,transfer:6,transmiss:[],transport:[3,4],transpos:6,tree:10,triangl:6,tril:6,tue:6,tup:6,turbul:[],two:6,txt:7,type:2,ucsd:7,unit:[3,4,6,10],univers:[4,7,10],unmodifi:3,until:5,updat:[5,6],upwel:6,url:7,usag:[2,6,7,10],use_banded_solv:3,user:[5,6],util:[],valid:2,valu:[5,6,7],vapor:10,vari:[4,7],variabl:[2,3,5,6],variat:7,vector:[6,7],veloc:7,verbos:[5,7],verison:6,vernal:7,version:6,vertic:[3,6],volumetr:6,walk:[],walk_process:10,water:[2,10],water_depth:[2,3,4,10],water_vapor:[],wavelength:10,wavenumb:10,weather:10,well:7,where:[3,7],whether:3,which:[2,3,5,7],within:10,wood:10,work:[3,7],www:7,yau:10,year:[5,7],zero:[5,6],zonal_mean_column:2,zonal_mean_surfac:[2,3]},titles:["climlab package","climlab.convection package","climlab.domain package","climlab.dynamics package","climlab.model package","climlab.process package","climlab.radiation package","climlab.solar package","climlab.surface package","climlab.tests package","climlab.utils package","Welcome to climlab-0.2.13&#8217;s documentation!","climlab-0.2.13","setup module"],titleterms:{albedo:8,aplusbt:6,axi:2,boltzmann:6,budyko_transport:3,climlab:[0,1,2,3,4,5,6,7,8,9,10,11,12],climtrad:6,cloud:6,column:4,constant:10,content:[0,1,2,3,4,5,6,7,8,9,10],convadj:1,convect:1,diagnost:5,diffus:3,document:11,domain:2,dynam:3,ebm:4,energy_budget:5,field:2,heat_capac:10,implicit:5,indic:11,insol:[6,7],legendr:10,model:4,modul:[0,1,2,3,4,5,6,7,8,9,10,13],nband:6,orbit:7,orbital_cycl:7,packag:[0,1,2,3,4,5,6,7,8,9,10],process:5,radiat:6,setup:13,solar:7,stommelbox:4,submodul:[1,2,3,4,5,6,7,8,10],subpackag:0,surfac:8,surface_radi:8,tabl:11,test:9,thermo:10,time_dependent_process:5,transmiss:6,turbul:8,util:10,walk:10,water_vapor:6,welcom:11}})