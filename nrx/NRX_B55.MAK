BORLANDPATH = "C:\program files\Borland\cbuilder5"

TASM = TASM32
TLIB = tlib
TLINK = ilink32
LIBPATH = $(BORLANDPATH)\LIB
INCLUDEPATH = $(BORLANDPATH)\INCLUDE

DIFF = sdiff
PRE =

CC = bcc32 -W- -v- -H- -3 -N -Og -Oi -Ov -f -I$(INCLUDEPATH)

.cpp.obj:
   $(CC) -c {$< }

everything:    	nrx_ex1.exe xcyclic.exe xgaulag.exe xpcshft.exe xshell.exe xairy.exe xdawson.exe xgauleg.exe xpearsn.exe xsimplx.exe xamebsa.exe xdbrent.exe xgaussj.exe xperiod.exe xsimpr.exe xamoeba.exe xddpoly.exe xgcf.exe xpiksr2.exe xsinft.exe xanneal.exe xdecchk.exe xgolden.exe xpiksrt.exe xsncndn.exe xarcode.exe xdf1dim.exe xgser.exe xplgndr.exe xsobseq.exe xavevar.exe xdfpmin.exe xhpsel.exe xpoidev.exe xsor.exe xbalanc.exe xdfridr.exe xhpsort.exe xpolcoe.exe xsort.exe xbandec.exe xdftint.exe xhqr.exe xpolcof.exe xsort2.exe xbanmul.exe xdlinmin.exe xhuffman.exe xpoldiv.exe xsort3.exe xbcucof.exe xeclass.exe xhunt.exe xpolin2.exe xspctrm.exe xbcuint.exe xeclazz.exe xhypgeo.exe xpolint.exe xspear.exe xbeschb.exe xei.exe xicrc.exe xpowell.exe xsphbes.exe xbessi.exe xeigsrt.exe xigray.exe xpredic.exe xbessi0.exe xelle.exe xindexx.exe xprobks.exe xsplie2.exe xbessi1.exe xellf.exe xirbit1.exe xpsdes.exe xsplin2.exe xbessik.exe xellpi.exe xirbit2.exe xpzextr.exe xspline.exe xbessj.exe xelmhes.exe xjacobi.exe xqgaus.exe xsplint.exe xbessj0.exe xerfcc.exe xjulday.exe xqrdcmp.exe xsprsax.exe xbessj1.exe xerff.exe xkendl1.exe xqromb.exe xsprsin.exe xbessjy.exe xerffc.exe xkendl2.exe xqromo.exe xsprspm.exe xbessk.exe xeulsum.exe xks2d1s.exe xqroot.exe xsprstm.exe xbessk0.exe xevlmem.exe xks2d2s.exe xqrsolv.exe xsprstp.exe xbessk1.exe xexpdev.exe xksone.exe xqrupdt.exe xsprstx.exe xbessy.exe xexpint.exe xkstwo.exe xqsimp.exe xstifbs.exe xbessy0.exe xf1dim.exe xlaguer.exe xqtrap.exe xstiff.exe xbessy1.exe xfactln.exe xlfit.exe xquad3d.exe xstoerm.exe xbeta.exe xfactrl.exe xlinbcg.exe xran.exe xsvbksb.exe xbetai.exe xfasper.exe xlinmin.exe xran4.exe xsvdcmp.exe xbico.exe xfgauss.exe xlocate.exe xrank.exe xsvdfit.exe xbnldev.exe xfit.exe xlubksb.exe xratint.exe xsvdvar.exe xbrent.exe xfitexy.exe xludcmp.exe xratlsq.exe xtoeplz.exe xbroydn.exe xfixrts.exe xmachar.exe xrc.exe xtptest.exe xbsstep.exe xfleg.exe xmedfit.exe xrd.exe xtqli.exe xcaldat.exe xflmoon.exe xmemcof.exe xrealft.exe xtrapzd.exe xchder.exe xfour1.exe xmgfas.exe xrf.exe xtred2.exe xchebev.exe xfourfs.exe xmglin.exe xrj.exe xtridag.exe xchebft.exe xfourn.exe xmidpnt.exe xrk4.exe xttest.exe xchebpc.exe xfpoly.exe xmiser.exe xrkdumb.exe xtutest.exe xchint.exe xfred2.exe xmmid.exe xrkqs.exe xtwofft.exe xcholsl.exe xfredin.exe xmnbrak.exe xrlft3.exe xvander.exe xchsone.exe xfrenel.exe xmnewt.exe xrofunc.exe xvegas.exe xchstwo.exe xfrprmn.exe xmoment.exe xrtbis.exe xvoltra.exe xcisi.exe xftest.exe xmppi.exe xrtflsp.exe xwt1.exe xcntab1.exe xgamdev.exe xmprove.exe xrtnewt.exe xwtn.exe xcntab2.exe xgammln.exe xmrqcof.exe xrtsafe.exe xzbrac.exe xconvlv.exe xgammp.exe xmrqmin.exe xrtsec.exe xzbrak.exe xcorrel.exe xgammq.exe xnewt.exe xrzextr.exe xzbrent.exe xcosft1.exe xgasdev.exe xodeint.exe xsavgol.exe xzrhqr.exe xcosft2.exe xgaucof.exe xorthog.exe xscrsho.exe xzriddr.exe xcovsrt.exe xgauher.exe xpade.exe xselect.exe xzroots.exe xcrank.exe xgaujac.exe xpccheb.exe xselip.exe 

nrx_lobj = addint.obj cntab2.obj frprmn.obj metrop.obj quad3d.obj snrm.obj airy.obj convlv.obj ftest.obj mgfas.obj quadct.obj sobseq.obj amebsa.obj copy.obj gamdev.obj mglin.obj quadmx.obj solvde.obj amoeba.obj correl.obj gammln.obj midexp.obj quadvl.obj sor.obj amotry.obj cosft1.obj gammp.obj midinf.obj ran0.obj sort1.obj amotsa.obj cosft2.obj gammq.obj midpnt.obj ran1.obj sort2.obj anneal.obj covsrt.obj gasdev.obj midsql.obj ran2.obj sort3.obj anorm2.obj crank.obj gaucof.obj midsqu.obj ran3.obj spctrm.obj arcmak.obj cyclic.obj gauher.obj miser.obj ran4.obj spear.obj arcode.obj daub4.obj gaujac.obj mmid.obj rank.obj sphbes.obj arcsum.obj dawson.obj gaulag.obj mnbrak.obj ranpt.obj sphfpt.obj asolve.obj dbrent.obj gauleg.obj mnewt.obj ratint.obj sphoot.obj atimes.obj ddpoly.obj gaussj.obj moment.obj ratlsq.obj splie2.obj avevar.obj decchk.obj gcf.obj mp2dfr.obj ratval.obj splin2.obj badluk.obj df1dim.obj golden.obj mpdiv.obj rc.obj spline.obj balanc.obj dfpmin.obj gser.obj mpinv.obj rd.obj splint.obj banbks.obj dfridr.obj hpsel.obj mpmul.obj realft.obj spread.obj bandec.obj dftcor.obj hpsort.obj mpops.obj rebin.obj sprsax.obj banmul.obj dftint.obj hqr.obj mppi.obj red.obj sprsin.obj bcucof.obj difeq.obj hufapp.obj mprove.obj relax.obj sprspm.obj bcuint.obj dlinmin.obj hufdec.obj mpsqrt.obj relax2.obj sprstm.obj beschb.obj eclass.obj hufenc.obj mrqcof.obj resid.obj sprstp.obj bessi.obj eclazz.obj hufmak.obj mrqmin.obj revcst.obj sprstx.obj bessi0.obj ei.obj hunt.obj newt.obj reverse.obj stifbs.obj bessi1.obj eigsrt.obj hypdrv.obj odeint.obj rf.obj stiff.obj bessik.obj elle.obj hypgeo.obj orthog.obj rj.obj stoerm.obj bessj.obj ellf.obj hypser.obj pade.obj rk4.obj svbksb.obj bessj0.obj ellpi.obj icrc.obj pccheb.obj rkck.obj svdcmp.obj bessj1.obj elmhes.obj icrc1.obj pcshft.obj rkdumb.obj svdfit.obj bessjy.obj erfcc.obj igray.obj pearsn.obj rkqs.obj svdvar.obj bessk.obj erff.obj indexx.obj period.obj rlft3.obj toeplz.obj bessk0.obj erffc.obj interp.obj piksr2.obj rofunc.obj tptest.obj bessk1.obj eulsum.obj irbit1.obj piksrt.obj rotate.obj tqli.obj bessy.obj evlmem.obj irbit2.obj pinvs.obj rsolv.obj trapzd.obj bessy0.obj expdev.obj jacobi1.obj plgndr.obj rstrct.obj tred2.obj bessy1.obj expint.obj jacobn.obj poidev.obj rtbis.obj tridag.obj beta.obj f1dim.obj julday.obj polcoe.obj rtflsp.obj trncst.obj betacf.obj factln.obj kendl1.obj polcof.obj rtnewt.obj trnspt.obj betai.obj factrl.obj kendl2.obj poldiv.obj rtsafe.obj ttest.obj bico.obj fasper.obj kermom.obj polin2.obj rtsec.obj tutest.obj bksub.obj fdjac.obj ks2d1s.obj polint.obj rzextr.obj twofft.obj bnldev.obj fgauss.obj ks2d2s.obj powell.obj savgol.obj vander.obj brent.obj ksone.obj predic.obj scrsho.obj vegas.obj broydn.obj fit.obj kstwo.obj probks.obj select.obj voltra.obj bsstep.obj fitexy.obj laguer.obj psdes.obj selip.obj wt1.obj caldat.obj fixrts.obj lfit.obj pwt.obj sfroid.obj wtn.obj chder.obj fleg.obj linbcg.obj pwtset.obj shell.obj wwghts.obj chebev.obj flmoon.obj linmin.obj pythag.obj shoot.obj zbrac.obj chebft.obj fmin.obj lnsrch.obj pzextr.obj shootf.obj zbrak.obj chebpc.obj four1.obj locate.obj qgaus.obj simp1.obj zbrent.obj chint.obj fourew.obj lop.obj qrdcmp.obj simp2.obj zrhqr.obj chixy.obj fourfs.obj lubksb.obj qromb.obj simp3.obj zriddr.obj choldc.obj fourn.obj ludcmp.obj qromo.obj simplx.obj zroots.obj cholsl.obj fpoly.obj machar.obj qroot.obj simpr.obj chsone.obj fred2.obj matadd.obj qrsolv.obj sinft.obj chstwo.obj fredex.obj matsub.obj qrupdt.obj slvsm2.obj cisi.obj fredin.obj medfit.obj qsimp.obj slvsml.obj cntab1.obj frenel.obj memcof.obj qtrap.obj sncndn.obj xsphfpt.obj nruti_nm.obj nrtyp_nm.obj newmat1.obj newmat2.obj newmat3.obj newmat4.obj newmat5.obj newmat6.obj newmat7.obj newmat8.obj newmatex.obj bandmat.obj submat.obj myexcept.obj

nrx.lib:       	$(nrx_lobj)
   $(TLIB) $@ /P32 /u $(nrx_lobj)

nrx_ex1_obj = nrx_ex1.obj cholesky.obj evalue.obj fft.obj hholder.obj jacobi.obj newfft.obj sort.obj svd.obj newmatrm.obj

nrx_ex1.exe:   	$(nrx_ex1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(nrx_ex1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcyclic_obj = xcyclic.obj

xcyclic.exe:   	$(xcyclic_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcyclic_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgaulag_obj = xgaulag.obj

xgaulag.exe:   	$(xgaulag_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgaulag_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpcshft_obj = xpcshft.obj

xpcshft.exe:   	$(xpcshft_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpcshft_obj),$@,, nrx.lib import32.lib cw32.lib
|

xshell_obj = xshell.obj

xshell.exe:    	$(xshell_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xshell_obj),$@,, nrx.lib import32.lib cw32.lib
|

xairy_obj = xairy.obj

xairy.exe:     	$(xairy_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xairy_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdawson_obj = xdawson.obj

xdawson.exe:   	$(xdawson_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdawson_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgauleg_obj = xgauleg.obj

xgauleg.exe:   	$(xgauleg_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgauleg_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpearsn_obj = xpearsn.obj

xpearsn.exe:   	$(xpearsn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpearsn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsimplx_obj = xsimplx.obj

xsimplx.exe:   	$(xsimplx_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsimplx_obj),$@,, nrx.lib import32.lib cw32.lib
|

xamebsa_obj = xamebsa.obj

xamebsa.exe:   	$(xamebsa_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xamebsa_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdbrent_obj = xdbrent.obj

xdbrent.exe:   	$(xdbrent_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdbrent_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgaussj_obj = xgaussj.obj

xgaussj.exe:   	$(xgaussj_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgaussj_obj),$@,, nrx.lib import32.lib cw32.lib
|

xperiod_obj = xperiod.obj

xperiod.exe:   	$(xperiod_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xperiod_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsimpr_obj = xsimpr.obj

xsimpr.exe:    	$(xsimpr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsimpr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xamoeba_obj = xamoeba.obj

xamoeba.exe:   	$(xamoeba_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xamoeba_obj),$@,, nrx.lib import32.lib cw32.lib
|

xddpoly_obj = xddpoly.obj

xddpoly.exe:   	$(xddpoly_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xddpoly_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgcf_obj = xgcf.obj

xgcf.exe:      	$(xgcf_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgcf_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpiksr2_obj = xpiksr2.obj

xpiksr2.exe:   	$(xpiksr2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpiksr2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsinft_obj = xsinft.obj

xsinft.exe:    	$(xsinft_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsinft_obj),$@,, nrx.lib import32.lib cw32.lib
|

xanneal_obj = xanneal.obj

xanneal.exe:   	$(xanneal_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xanneal_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdecchk_obj = xdecchk.obj

xdecchk.exe:   	$(xdecchk_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdecchk_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgolden_obj = xgolden.obj

xgolden.exe:   	$(xgolden_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgolden_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpiksrt_obj = xpiksrt.obj

xpiksrt.exe:   	$(xpiksrt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpiksrt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsncndn_obj = xsncndn.obj

xsncndn.exe:   	$(xsncndn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsncndn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xarcode_obj = xarcode.obj

xarcode.exe:   	$(xarcode_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xarcode_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdf1dim_obj = xdf1dim.obj

xdf1dim.exe:   	$(xdf1dim_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdf1dim_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgser_obj = xgser.obj

xgser.exe:     	$(xgser_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgser_obj),$@,, nrx.lib import32.lib cw32.lib
|

xplgndr_obj = xplgndr.obj

xplgndr.exe:   	$(xplgndr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xplgndr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsobseq_obj = xsobseq.obj

xsobseq.exe:   	$(xsobseq_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsobseq_obj),$@,, nrx.lib import32.lib cw32.lib
|

xavevar_obj = xavevar.obj

xavevar.exe:   	$(xavevar_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xavevar_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdfpmin_obj = xdfpmin.obj

xdfpmin.exe:   	$(xdfpmin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdfpmin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhpsel_obj = xhpsel.obj

xhpsel.exe:    	$(xhpsel_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhpsel_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpoidev_obj = xpoidev.obj

xpoidev.exe:   	$(xpoidev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpoidev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsor_obj = xsor.obj

xsor.exe:      	$(xsor_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsor_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbalanc_obj = xbalanc.obj

xbalanc.exe:   	$(xbalanc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbalanc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdfridr_obj = xdfridr.obj

xdfridr.exe:   	$(xdfridr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdfridr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhpsort_obj = xhpsort.obj

xhpsort.exe:   	$(xhpsort_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhpsort_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpolcoe_obj = xpolcoe.obj

xpolcoe.exe:   	$(xpolcoe_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpolcoe_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsort_obj = xsort.obj

xsort.exe:     	$(xsort_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsort_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbandec_obj = xbandec.obj

xbandec.exe:   	$(xbandec_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbandec_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdftint_obj = xdftint.obj

xdftint.exe:   	$(xdftint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdftint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhqr_obj = xhqr.obj

xhqr.exe:      	$(xhqr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhqr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpolcof_obj = xpolcof.obj

xpolcof.exe:   	$(xpolcof_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpolcof_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsort2_obj = xsort2.obj

xsort2.exe:    	$(xsort2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsort2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbanmul_obj = xbanmul.obj

xbanmul.exe:   	$(xbanmul_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbanmul_obj),$@,, nrx.lib import32.lib cw32.lib
|

xdlinmin_obj = xdlinmin.obj

xdlinmin.exe:  	$(xdlinmin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xdlinmin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhuffman_obj = xhuffman.obj

xhuffman.exe:  	$(xhuffman_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhuffman_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpoldiv_obj = xpoldiv.obj

xpoldiv.exe:   	$(xpoldiv_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpoldiv_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsort3_obj = xsort3.obj

xsort3.exe:    	$(xsort3_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsort3_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbcucof_obj = xbcucof.obj

xbcucof.exe:   	$(xbcucof_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbcucof_obj),$@,, nrx.lib import32.lib cw32.lib
|

xeclass_obj = xeclass.obj

xeclass.exe:   	$(xeclass_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xeclass_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhunt_obj = xhunt.obj

xhunt.exe:     	$(xhunt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhunt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpolin2_obj = xpolin2.obj

xpolin2.exe:   	$(xpolin2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpolin2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xspctrm_obj = xspctrm.obj

xspctrm.exe:   	$(xspctrm_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xspctrm_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbcuint_obj = xbcuint.obj

xbcuint.exe:   	$(xbcuint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbcuint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xeclazz_obj = xeclazz.obj

xeclazz.exe:   	$(xeclazz_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xeclazz_obj),$@,, nrx.lib import32.lib cw32.lib
|

xhypgeo_obj = xhypgeo.obj

xhypgeo.exe:   	$(xhypgeo_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xhypgeo_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpolint_obj = xpolint.obj

xpolint.exe:   	$(xpolint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpolint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xspear_obj = xspear.obj

xspear.exe:    	$(xspear_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xspear_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbeschb_obj = xbeschb.obj

xbeschb.exe:   	$(xbeschb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbeschb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xei_obj = xei.obj

xei.exe:       	$(xei_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xei_obj),$@,, nrx.lib import32.lib cw32.lib
|

xicrc_obj = xicrc.obj

xicrc.exe:     	$(xicrc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xicrc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpowell_obj = xpowell.obj

xpowell.exe:   	$(xpowell_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpowell_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsphbes_obj = xsphbes.obj

xsphbes.exe:   	$(xsphbes_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsphbes_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessi_obj = xbessi.obj

xbessi.exe:    	$(xbessi_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessi_obj),$@,, nrx.lib import32.lib cw32.lib
|

xeigsrt_obj = xeigsrt.obj

xeigsrt.exe:   	$(xeigsrt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xeigsrt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xigray_obj = xigray.obj

xigray.exe:    	$(xigray_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xigray_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpredic_obj = xpredic.obj

xpredic.exe:   	$(xpredic_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpredic_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessi0_obj = xbessi0.obj

xbessi0.exe:   	$(xbessi0_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessi0_obj),$@,, nrx.lib import32.lib cw32.lib
|

xelle_obj = xelle.obj

xelle.exe:     	$(xelle_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xelle_obj),$@,, nrx.lib import32.lib cw32.lib
|

xindexx_obj = xindexx.obj

xindexx.exe:   	$(xindexx_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xindexx_obj),$@,, nrx.lib import32.lib cw32.lib
|

xprobks_obj = xprobks.obj

xprobks.exe:   	$(xprobks_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xprobks_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsplie2_obj = xsplie2.obj

xsplie2.exe:   	$(xsplie2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsplie2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessi1_obj = xbessi1.obj

xbessi1.exe:   	$(xbessi1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessi1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xellf_obj = xellf.obj

xellf.exe:     	$(xellf_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xellf_obj),$@,, nrx.lib import32.lib cw32.lib
|

xirbit1_obj = xirbit1.obj

xirbit1.exe:   	$(xirbit1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xirbit1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpsdes_obj = xpsdes.obj

xpsdes.exe:    	$(xpsdes_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpsdes_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsplin2_obj = xsplin2.obj

xsplin2.exe:   	$(xsplin2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsplin2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessik_obj = xbessik.obj

xbessik.exe:   	$(xbessik_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessik_obj),$@,, nrx.lib import32.lib cw32.lib
|

xellpi_obj = xellpi.obj

xellpi.exe:    	$(xellpi_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xellpi_obj),$@,, nrx.lib import32.lib cw32.lib
|

xirbit2_obj = xirbit2.obj

xirbit2.exe:   	$(xirbit2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xirbit2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpzextr_obj = xpzextr.obj

xpzextr.exe:   	$(xpzextr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpzextr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xspline_obj = xspline.obj

xspline.exe:   	$(xspline_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xspline_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessj_obj = xbessj.obj

xbessj.exe:    	$(xbessj_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessj_obj),$@,, nrx.lib import32.lib cw32.lib
|

xelmhes_obj = xelmhes.obj

xelmhes.exe:   	$(xelmhes_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xelmhes_obj),$@,, nrx.lib import32.lib cw32.lib
|

xjacobi_obj = xjacobi.obj

xjacobi.exe:   	$(xjacobi_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xjacobi_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqgaus_obj = xqgaus.obj

xqgaus.exe:    	$(xqgaus_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqgaus_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsplint_obj = xsplint.obj

xsplint.exe:   	$(xsplint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsplint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessj0_obj = xbessj0.obj

xbessj0.exe:   	$(xbessj0_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessj0_obj),$@,, nrx.lib import32.lib cw32.lib
|

xerfcc_obj = xerfcc.obj

xerfcc.exe:    	$(xerfcc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xerfcc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xjulday_obj = xjulday.obj

xjulday.exe:   	$(xjulday_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xjulday_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqrdcmp_obj = xqrdcmp.obj

xqrdcmp.exe:   	$(xqrdcmp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqrdcmp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprsax_obj = xsprsax.obj

xsprsax.exe:   	$(xsprsax_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprsax_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessj1_obj = xbessj1.obj

xbessj1.exe:   	$(xbessj1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessj1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xerff_obj = xerff.obj

xerff.exe:     	$(xerff_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xerff_obj),$@,, nrx.lib import32.lib cw32.lib
|

xkendl1_obj = xkendl1.obj

xkendl1.exe:   	$(xkendl1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xkendl1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqromb_obj = xqromb.obj

xqromb.exe:    	$(xqromb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqromb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprsin_obj = xsprsin.obj

xsprsin.exe:   	$(xsprsin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprsin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessjy_obj = xbessjy.obj

xbessjy.exe:   	$(xbessjy_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessjy_obj),$@,, nrx.lib import32.lib cw32.lib
|

xerffc_obj = xerffc.obj

xerffc.exe:    	$(xerffc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xerffc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xkendl2_obj = xkendl2.obj

xkendl2.exe:   	$(xkendl2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xkendl2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqromo_obj = xqromo.obj

xqromo.exe:    	$(xqromo_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqromo_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprspm_obj = xsprspm.obj

xsprspm.exe:   	$(xsprspm_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprspm_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessk_obj = xbessk.obj

xbessk.exe:    	$(xbessk_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessk_obj),$@,, nrx.lib import32.lib cw32.lib
|

xeulsum_obj = xeulsum.obj

xeulsum.exe:   	$(xeulsum_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xeulsum_obj),$@,, nrx.lib import32.lib cw32.lib
|

xks2d1s_obj = xks2d1s.obj

xks2d1s.exe:   	$(xks2d1s_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xks2d1s_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqroot_obj = xqroot.obj

xqroot.exe:    	$(xqroot_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqroot_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprstm_obj = xsprstm.obj

xsprstm.exe:   	$(xsprstm_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprstm_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessk0_obj = xbessk0.obj

xbessk0.exe:   	$(xbessk0_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessk0_obj),$@,, nrx.lib import32.lib cw32.lib
|

xevlmem_obj = xevlmem.obj

xevlmem.exe:   	$(xevlmem_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xevlmem_obj),$@,, nrx.lib import32.lib cw32.lib
|

xks2d2s_obj = xks2d2s.obj

xks2d2s.exe:   	$(xks2d2s_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xks2d2s_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqrsolv_obj = xqrsolv.obj

xqrsolv.exe:   	$(xqrsolv_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqrsolv_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprstp_obj = xsprstp.obj

xsprstp.exe:   	$(xsprstp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprstp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessk1_obj = xbessk1.obj

xbessk1.exe:   	$(xbessk1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessk1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xexpdev_obj = xexpdev.obj

xexpdev.exe:   	$(xexpdev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xexpdev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xksone_obj = xksone.obj

xksone.exe:    	$(xksone_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xksone_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqrupdt_obj = xqrupdt.obj

xqrupdt.exe:   	$(xqrupdt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqrupdt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsprstx_obj = xsprstx.obj

xsprstx.exe:   	$(xsprstx_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsprstx_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessy_obj = xbessy.obj

xbessy.exe:    	$(xbessy_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessy_obj),$@,, nrx.lib import32.lib cw32.lib
|

xexpint_obj = xexpint.obj

xexpint.exe:   	$(xexpint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xexpint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xkstwo_obj = xkstwo.obj

xkstwo.exe:    	$(xkstwo_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xkstwo_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqsimp_obj = xqsimp.obj

xqsimp.exe:    	$(xqsimp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqsimp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xstifbs_obj = xstifbs.obj

xstifbs.exe:   	$(xstifbs_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xstifbs_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessy0_obj = xbessy0.obj

xbessy0.exe:   	$(xbessy0_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessy0_obj),$@,, nrx.lib import32.lib cw32.lib
|

xf1dim_obj = xf1dim.obj

xf1dim.exe:    	$(xf1dim_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xf1dim_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlaguer_obj = xlaguer.obj

xlaguer.exe:   	$(xlaguer_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlaguer_obj),$@,, nrx.lib import32.lib cw32.lib
|

xqtrap_obj = xqtrap.obj

xqtrap.exe:    	$(xqtrap_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xqtrap_obj),$@,, nrx.lib import32.lib cw32.lib
|

xstiff_obj = xstiff.obj

xstiff.exe:    	$(xstiff_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xstiff_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbessy1_obj = xbessy1.obj

xbessy1.exe:   	$(xbessy1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbessy1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfactln_obj = xfactln.obj

xfactln.exe:   	$(xfactln_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfactln_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlfit_obj = xlfit.obj

xlfit.exe:     	$(xlfit_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlfit_obj),$@,, nrx.lib import32.lib cw32.lib
|

xquad3d_obj = xquad3d.obj

xquad3d.exe:   	$(xquad3d_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xquad3d_obj),$@,, nrx.lib import32.lib cw32.lib
|

xstoerm_obj = xstoerm.obj

xstoerm.exe:   	$(xstoerm_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xstoerm_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbeta_obj = xbeta.obj

xbeta.exe:     	$(xbeta_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbeta_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfactrl_obj = xfactrl.obj

xfactrl.exe:   	$(xfactrl_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfactrl_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlinbcg_obj = xlinbcg.obj

xlinbcg.exe:   	$(xlinbcg_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlinbcg_obj),$@,, nrx.lib import32.lib cw32.lib
|

xran_obj = xran.obj

xran.exe:      	$(xran_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xran_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsvbksb_obj = xsvbksb.obj

xsvbksb.exe:   	$(xsvbksb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsvbksb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbetai_obj = xbetai.obj

xbetai.exe:    	$(xbetai_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbetai_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfasper_obj = xfasper.obj

xfasper.exe:   	$(xfasper_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfasper_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlinmin_obj = xlinmin.obj

xlinmin.exe:   	$(xlinmin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlinmin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xran4_obj = xran4.obj

xran4.exe:     	$(xran4_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xran4_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsvdcmp_obj = xsvdcmp.obj

xsvdcmp.exe:   	$(xsvdcmp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsvdcmp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbico_obj = xbico.obj

xbico.exe:     	$(xbico_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbico_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfgauss_obj = xfgauss.obj

xfgauss.exe:   	$(xfgauss_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfgauss_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlocate_obj = xlocate.obj

xlocate.exe:   	$(xlocate_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlocate_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrank_obj = xrank.obj

xrank.exe:     	$(xrank_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrank_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsvdfit_obj = xsvdfit.obj

xsvdfit.exe:   	$(xsvdfit_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsvdfit_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbnldev_obj = xbnldev.obj

xbnldev.exe:   	$(xbnldev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbnldev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfit_obj = xfit.obj

xfit.exe:      	$(xfit_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfit_obj),$@,, nrx.lib import32.lib cw32.lib
|

xlubksb_obj = xlubksb.obj

xlubksb.exe:   	$(xlubksb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xlubksb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xratint_obj = xratint.obj

xratint.exe:   	$(xratint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xratint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsvdvar_obj = xsvdvar.obj

xsvdvar.exe:   	$(xsvdvar_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsvdvar_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbrent_obj = xbrent.obj

xbrent.exe:    	$(xbrent_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbrent_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfitexy_obj = xfitexy.obj

xfitexy.exe:   	$(xfitexy_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfitexy_obj),$@,, nrx.lib import32.lib cw32.lib
|

xludcmp_obj = xludcmp.obj

xludcmp.exe:   	$(xludcmp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xludcmp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xratlsq_obj = xratlsq.obj

xratlsq.exe:   	$(xratlsq_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xratlsq_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtoeplz_obj = xtoeplz.obj

xtoeplz.exe:   	$(xtoeplz_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtoeplz_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbroydn_obj = xbroydn.obj

xbroydn.exe:   	$(xbroydn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbroydn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfixrts_obj = xfixrts.obj

xfixrts.exe:   	$(xfixrts_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfixrts_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmachar_obj = xmachar.obj

xmachar.exe:   	$(xmachar_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmachar_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrc_obj = xrc.obj

xrc.exe:       	$(xrc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtptest_obj = xtptest.obj

xtptest.exe:   	$(xtptest_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtptest_obj),$@,, nrx.lib import32.lib cw32.lib
|

xbsstep_obj = xbsstep.obj

xbsstep.exe:   	$(xbsstep_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xbsstep_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfleg_obj = xfleg.obj

xfleg.exe:     	$(xfleg_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfleg_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmedfit_obj = xmedfit.obj

xmedfit.exe:   	$(xmedfit_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmedfit_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrd_obj = xrd.obj

xrd.exe:       	$(xrd_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrd_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtqli_obj = xtqli.obj

xtqli.exe:     	$(xtqli_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtqli_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcaldat_obj = xcaldat.obj

xcaldat.exe:   	$(xcaldat_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcaldat_obj),$@,, nrx.lib import32.lib cw32.lib
|

xflmoon_obj = xflmoon.obj

xflmoon.exe:   	$(xflmoon_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xflmoon_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmemcof_obj = xmemcof.obj

xmemcof.exe:   	$(xmemcof_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmemcof_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrealft_obj = xrealft.obj

xrealft.exe:   	$(xrealft_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrealft_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtrapzd_obj = xtrapzd.obj

xtrapzd.exe:   	$(xtrapzd_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtrapzd_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchder_obj = xchder.obj

xchder.exe:    	$(xchder_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchder_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfour1_obj = xfour1.obj

xfour1.exe:    	$(xfour1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfour1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmgfas_obj = xmgfas.obj

xmgfas.exe:    	$(xmgfas_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmgfas_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrf_obj = xrf.obj

xrf.exe:       	$(xrf_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrf_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtred2_obj = xtred2.obj

xtred2.exe:    	$(xtred2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtred2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchebev_obj = xchebev.obj

xchebev.exe:   	$(xchebev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchebev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfourfs_obj = xfourfs.obj

xfourfs.exe:   	$(xfourfs_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfourfs_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmglin_obj = xmglin.obj

xmglin.exe:    	$(xmglin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmglin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrj_obj = xrj.obj

xrj.exe:       	$(xrj_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrj_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtridag_obj = xtridag.obj

xtridag.exe:   	$(xtridag_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtridag_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchebft_obj = xchebft.obj

xchebft.exe:   	$(xchebft_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchebft_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfourn_obj = xfourn.obj

xfourn.exe:    	$(xfourn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfourn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmidpnt_obj = xmidpnt.obj

xmidpnt.exe:   	$(xmidpnt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmidpnt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrk4_obj = xrk4.obj

xrk4.exe:      	$(xrk4_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrk4_obj),$@,, nrx.lib import32.lib cw32.lib
|

xttest_obj = xttest.obj

xttest.exe:    	$(xttest_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xttest_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchebpc_obj = xchebpc.obj

xchebpc.exe:   	$(xchebpc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchebpc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfpoly_obj = xfpoly.obj

xfpoly.exe:    	$(xfpoly_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfpoly_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmiser_obj = xmiser.obj

xmiser.exe:    	$(xmiser_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmiser_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrkdumb_obj = xrkdumb.obj

xrkdumb.exe:   	$(xrkdumb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrkdumb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtutest_obj = xtutest.obj

xtutest.exe:   	$(xtutest_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtutest_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchint_obj = xchint.obj

xchint.exe:    	$(xchint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfred2_obj = xfred2.obj

xfred2.exe:    	$(xfred2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfred2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmmid_obj = xmmid.obj

xmmid.exe:     	$(xmmid_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmmid_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrkqs_obj = xrkqs.obj

xrkqs.exe:     	$(xrkqs_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrkqs_obj),$@,, nrx.lib import32.lib cw32.lib
|

xtwofft_obj = xtwofft.obj

xtwofft.exe:   	$(xtwofft_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xtwofft_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcholsl_obj = xcholsl.obj

xcholsl.exe:   	$(xcholsl_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcholsl_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfredin_obj = xfredin.obj

xfredin.exe:   	$(xfredin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfredin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmnbrak_obj = xmnbrak.obj

xmnbrak.exe:   	$(xmnbrak_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmnbrak_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrlft3_obj = xrlft3.obj

xrlft3.exe:    	$(xrlft3_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrlft3_obj),$@,, nrx.lib import32.lib cw32.lib
|

xvander_obj = xvander.obj

xvander.exe:   	$(xvander_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xvander_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchsone_obj = xchsone.obj

xchsone.exe:   	$(xchsone_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchsone_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfrenel_obj = xfrenel.obj

xfrenel.exe:   	$(xfrenel_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfrenel_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmnewt_obj = xmnewt.obj

xmnewt.exe:    	$(xmnewt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmnewt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrofunc_obj = xrofunc.obj

xrofunc.exe:   	$(xrofunc_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrofunc_obj),$@,, nrx.lib import32.lib cw32.lib
|

xvegas_obj = xvegas.obj

xvegas.exe:    	$(xvegas_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xvegas_obj),$@,, nrx.lib import32.lib cw32.lib
|

xchstwo_obj = xchstwo.obj

xchstwo.exe:   	$(xchstwo_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xchstwo_obj),$@,, nrx.lib import32.lib cw32.lib
|

xfrprmn_obj = xfrprmn.obj

xfrprmn.exe:   	$(xfrprmn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xfrprmn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmoment_obj = xmoment.obj

xmoment.exe:   	$(xmoment_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmoment_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrtbis_obj = xrtbis.obj

xrtbis.exe:    	$(xrtbis_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrtbis_obj),$@,, nrx.lib import32.lib cw32.lib
|

xvoltra_obj = xvoltra.obj

xvoltra.exe:   	$(xvoltra_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xvoltra_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcisi_obj = xcisi.obj

xcisi.exe:     	$(xcisi_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcisi_obj),$@,, nrx.lib import32.lib cw32.lib
|

xftest_obj = xftest.obj

xftest.exe:    	$(xftest_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xftest_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmppi_obj = xmppi.obj

xmppi.exe:     	$(xmppi_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmppi_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrtflsp_obj = xrtflsp.obj

xrtflsp.exe:   	$(xrtflsp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrtflsp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xwt1_obj = xwt1.obj

xwt1.exe:      	$(xwt1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xwt1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcntab1_obj = xcntab1.obj

xcntab1.exe:   	$(xcntab1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcntab1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgamdev_obj = xgamdev.obj

xgamdev.exe:   	$(xgamdev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgamdev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmprove_obj = xmprove.obj

xmprove.exe:   	$(xmprove_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmprove_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrtnewt_obj = xrtnewt.obj

xrtnewt.exe:   	$(xrtnewt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrtnewt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xwtn_obj = xwtn.obj

xwtn.exe:      	$(xwtn_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xwtn_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcntab2_obj = xcntab2.obj

xcntab2.exe:   	$(xcntab2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcntab2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgammln_obj = xgammln.obj

xgammln.exe:   	$(xgammln_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgammln_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmrqcof_obj = xmrqcof.obj

xmrqcof.exe:   	$(xmrqcof_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmrqcof_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrtsafe_obj = xrtsafe.obj

xrtsafe.exe:   	$(xrtsafe_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrtsafe_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzbrac_obj = xzbrac.obj

xzbrac.exe:    	$(xzbrac_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzbrac_obj),$@,, nrx.lib import32.lib cw32.lib
|

xconvlv_obj = xconvlv.obj

xconvlv.exe:   	$(xconvlv_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xconvlv_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgammp_obj = xgammp.obj

xgammp.exe:    	$(xgammp_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgammp_obj),$@,, nrx.lib import32.lib cw32.lib
|

xmrqmin_obj = xmrqmin.obj

xmrqmin.exe:   	$(xmrqmin_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xmrqmin_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrtsec_obj = xrtsec.obj

xrtsec.exe:    	$(xrtsec_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrtsec_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzbrak_obj = xzbrak.obj

xzbrak.exe:    	$(xzbrak_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzbrak_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcorrel_obj = xcorrel.obj

xcorrel.exe:   	$(xcorrel_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcorrel_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgammq_obj = xgammq.obj

xgammq.exe:    	$(xgammq_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgammq_obj),$@,, nrx.lib import32.lib cw32.lib
|

xnewt_obj = xnewt.obj

xnewt.exe:     	$(xnewt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xnewt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xrzextr_obj = xrzextr.obj

xrzextr.exe:   	$(xrzextr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xrzextr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzbrent_obj = xzbrent.obj

xzbrent.exe:   	$(xzbrent_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzbrent_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcosft1_obj = xcosft1.obj

xcosft1.exe:   	$(xcosft1_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcosft1_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgasdev_obj = xgasdev.obj

xgasdev.exe:   	$(xgasdev_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgasdev_obj),$@,, nrx.lib import32.lib cw32.lib
|

xodeint_obj = xodeint.obj

xodeint.exe:   	$(xodeint_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xodeint_obj),$@,, nrx.lib import32.lib cw32.lib
|

xsavgol_obj = xsavgol.obj

xsavgol.exe:   	$(xsavgol_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xsavgol_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzrhqr_obj = xzrhqr.obj

xzrhqr.exe:    	$(xzrhqr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzrhqr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcosft2_obj = xcosft2.obj

xcosft2.exe:   	$(xcosft2_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcosft2_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgaucof_obj = xgaucof.obj

xgaucof.exe:   	$(xgaucof_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgaucof_obj),$@,, nrx.lib import32.lib cw32.lib
|

xorthog_obj = xorthog.obj

xorthog.exe:   	$(xorthog_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xorthog_obj),$@,, nrx.lib import32.lib cw32.lib
|

xscrsho_obj = xscrsho.obj

xscrsho.exe:   	$(xscrsho_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xscrsho_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzriddr_obj = xzriddr.obj

xzriddr.exe:   	$(xzriddr_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzriddr_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcovsrt_obj = xcovsrt.obj

xcovsrt.exe:   	$(xcovsrt_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcovsrt_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgauher_obj = xgauher.obj

xgauher.exe:   	$(xgauher_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgauher_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpade_obj = xpade.obj

xpade.exe:     	$(xpade_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpade_obj),$@,, nrx.lib import32.lib cw32.lib
|

xselect_obj = xselect.obj

xselect.exe:   	$(xselect_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xselect_obj),$@,, nrx.lib import32.lib cw32.lib
|

xzroots_obj = xzroots.obj

xzroots.exe:   	$(xzroots_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xzroots_obj),$@,, nrx.lib import32.lib cw32.lib
|

xcrank_obj = xcrank.obj

xcrank.exe:    	$(xcrank_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xcrank_obj),$@,, nrx.lib import32.lib cw32.lib
|

xgaujac_obj = xgaujac.obj

xgaujac.exe:   	$(xgaujac_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xgaujac_obj),$@,, nrx.lib import32.lib cw32.lib
|

xpccheb_obj = xpccheb.obj

xpccheb.exe:   	$(xpccheb_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xpccheb_obj),$@,, nrx.lib import32.lib cw32.lib
|

xselip_obj = xselip.obj

xselip.exe:    	$(xselip_obj) nrx.lib
   $(TLINK) /x/L$(LIBPATH)/Gn -Tpe -ap -c @&&|
c0x32.obj $(xselip_obj),$@,, nrx.lib import32.lib cw32.lib
|

addint.obj:    	addint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cntab2.obj:    	cntab2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

frprmn.obj:    	frprmn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

metrop.obj:    	metrop.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quad3d.obj:    	quad3d.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

snrm.obj:      	snrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

airy.obj:      	airy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

convlv.obj:    	convlv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ftest.obj:     	ftest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mgfas.obj:     	mgfas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadct.obj:    	quadct.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sobseq.obj:    	sobseq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amebsa.obj:    	amebsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

copy.obj:      	copy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gamdev.obj:    	gamdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mglin.obj:     	mglin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadmx.obj:    	quadmx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

solvde.obj:    	solvde.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amoeba.obj:    	amoeba.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

correl.obj:    	correl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammln.obj:    	gammln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midexp.obj:    	midexp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadvl.obj:    	quadvl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sor.obj:       	sor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amotry.obj:    	amotry.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cosft1.obj:    	cosft1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammp.obj:     	gammp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midinf.obj:    	midinf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran0.obj:      	ran0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort1.obj:     	sort1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amotsa.obj:    	amotsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cosft2.obj:    	cosft2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammq.obj:     	gammq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midpnt.obj:    	midpnt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran1.obj:      	ran1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort2.obj:     	sort2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

anneal.obj:    	anneal.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

covsrt.obj:    	covsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gasdev.obj:    	gasdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midsql.obj:    	midsql.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran2.obj:      	ran2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort3.obj:     	sort3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

anorm2.obj:    	anorm2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

crank.obj:     	crank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaucof.obj:    	gaucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midsqu.obj:    	midsqu.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran3.obj:      	ran3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spctrm.obj:    	spctrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcmak.obj:    	arcmak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cyclic.obj:    	cyclic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gauher.obj:    	gauher.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

miser.obj:     	miser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran4.obj:      	ran4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spear.obj:     	spear.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcode.obj:    	arcode.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

daub4.obj:     	daub4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaujac.obj:    	gaujac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mmid.obj:      	mmid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rank.obj:      	rank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphbes.obj:    	sphbes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcsum.obj:    	arcsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dawson.obj:    	dawson.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaulag.obj:    	gaulag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mnbrak.obj:    	mnbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ranpt.obj:     	ranpt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphfpt.obj:    	sphfpt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

asolve.obj:    	asolve.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dbrent.obj:    	dbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gauleg.obj:    	gauleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mnewt.obj:     	mnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratint.obj:    	ratint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphoot.obj:    	sphoot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

atimes.obj:    	atimes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ddpoly.obj:    	ddpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaussj.obj:    	gaussj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

moment.obj:    	moment.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratlsq.obj:    	ratlsq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splie2.obj:    	splie2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

avevar.obj:    	avevar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

decchk.obj:    	decchk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gcf.obj:       	gcf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mp2dfr.obj:    	mp2dfr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratval.obj:    	ratval.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splin2.obj:    	splin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

badluk.obj:    	badluk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

df1dim.obj:    	df1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

golden.obj:    	golden.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpdiv.obj:     	mpdiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rc.obj:        	rc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spline.obj:    	spline.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

balanc.obj:    	balanc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dfpmin.obj:    	dfpmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gser.obj:      	gser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpinv.obj:     	mpinv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rd.obj:        	rd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splint.obj:    	splint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

banbks.obj:    	banbks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dfridr.obj:    	dfridr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hpsel.obj:     	hpsel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpmul.obj:     	mpmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

realft.obj:    	realft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spread.obj:    	spread.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bandec.obj:    	bandec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dftcor.obj:    	dftcor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hpsort.obj:    	hpsort.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpops.obj:     	mpops.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rebin.obj:     	rebin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprsax.obj:    	sprsax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

banmul.obj:    	banmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dftint.obj:    	dftint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hqr.obj:       	hqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mppi.obj:      	mppi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

red.obj:       	red.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprsin.obj:    	sprsin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bcucof.obj:    	bcucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

difeq.obj:     	difeq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufapp.obj:    	hufapp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mprove.obj:    	mprove.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

relax.obj:     	relax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprspm.obj:    	sprspm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bcuint.obj:    	bcuint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dlinmin.obj:   	dlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufdec.obj:    	hufdec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpsqrt.obj:    	mpsqrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

relax2.obj:    	relax2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstm.obj:    	sprstm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

beschb.obj:    	beschb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eclass.obj:    	eclass.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufenc.obj:    	hufenc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mrqcof.obj:    	mrqcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

resid.obj:     	resid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstp.obj:    	sprstp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi.obj:     	bessi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eclazz.obj:    	eclazz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufmak.obj:    	hufmak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mrqmin.obj:    	mrqmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

revcst.obj:    	revcst.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstx.obj:    	sprstx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi0.obj:    	bessi0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ei.obj:        	ei.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hunt.obj:      	hunt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

newt.obj:      	newt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

reverse.obj:   	reverse.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stifbs.obj:    	stifbs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi1.obj:    	bessi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eigsrt.obj:    	eigsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypdrv.obj:    	hypdrv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

odeint.obj:    	odeint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rf.obj:        	rf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stiff.obj:     	stiff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessik.obj:    	bessik.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

elle.obj:      	elle.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypgeo.obj:    	hypgeo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

orthog.obj:    	orthog.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rj.obj:        	rj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stoerm.obj:    	stoerm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj.obj:     	bessj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ellf.obj:      	ellf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypser.obj:    	hypser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pade.obj:      	pade.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rk4.obj:       	rk4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svbksb.obj:    	svbksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj0.obj:    	bessj0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ellpi.obj:     	ellpi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

icrc.obj:      	icrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pccheb.obj:    	pccheb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkck.obj:      	rkck.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdcmp.obj:    	svdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj1.obj:    	bessj1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

elmhes.obj:    	elmhes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

icrc1.obj:     	icrc1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pcshft.obj:    	pcshft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkdumb.obj:    	rkdumb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdfit.obj:    	svdfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessjy.obj:    	bessjy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erfcc.obj:     	erfcc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

igray.obj:     	igray.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pearsn.obj:    	pearsn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkqs.obj:      	rkqs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdvar.obj:    	svdvar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk.obj:     	bessk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erff.obj:      	erff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

indexx.obj:    	indexx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

period.obj:    	period.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rlft3.obj:     	rlft3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

toeplz.obj:    	toeplz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk0.obj:    	bessk0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erffc.obj:     	erffc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

interp.obj:    	interp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

piksr2.obj:    	piksr2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rofunc.obj:    	rofunc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tptest.obj:    	tptest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk1.obj:    	bessk1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eulsum.obj:    	eulsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

irbit1.obj:    	irbit1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

piksrt.obj:    	piksrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rotate.obj:    	rotate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tqli.obj:      	tqli.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy.obj:     	bessy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

evlmem.obj:    	evlmem.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

irbit2.obj:    	irbit2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pinvs.obj:     	pinvs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rsolv.obj:     	rsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trapzd.obj:    	trapzd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy0.obj:    	bessy0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

expdev.obj:    	expdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

jacobi1.obj:   	jacobi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

plgndr.obj:    	plgndr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rstrct.obj:    	rstrct.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tred2.obj:     	tred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy1.obj:    	bessy1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

expint.obj:    	expint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

jacobn.obj:    	jacobn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

poidev.obj:    	poidev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtbis.obj:     	rtbis.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tridag.obj:    	tridag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

beta.obj:      	beta.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

f1dim.obj:     	f1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

julday.obj:    	julday.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polcoe.obj:    	polcoe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtflsp.obj:    	rtflsp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trncst.obj:    	trncst.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

betacf.obj:    	betacf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

factln.obj:    	factln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kendl1.obj:    	kendl1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polcof.obj:    	polcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtnewt.obj:    	rtnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trnspt.obj:    	trnspt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

betai.obj:     	betai.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

factrl.obj:    	factrl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kendl2.obj:    	kendl2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

poldiv.obj:    	poldiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtsafe.obj:    	rtsafe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ttest.obj:     	ttest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bico.obj:      	bico.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fasper.obj:    	fasper.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kermom.obj:    	kermom.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polin2.obj:    	polin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtsec.obj:     	rtsec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tutest.obj:    	tutest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bksub.obj:     	bksub.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fdjac.obj:     	fdjac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ks2d1s.obj:    	ks2d1s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polint.obj:    	polint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rzextr.obj:    	rzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

twofft.obj:    	twofft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bnldev.obj:    	bnldev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fgauss.obj:    	fgauss.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ks2d2s.obj:    	ks2d2s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

powell.obj:    	powell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

savgol.obj:    	savgol.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

vander.obj:    	vander.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

brent.obj:     	brent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ksone.obj:     	ksone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

predic.obj:    	predic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

scrsho.obj:    	scrsho.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

vegas.obj:     	vegas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

broydn.obj:    	broydn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fit.obj:       	fit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kstwo.obj:     	kstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

probks.obj:    	probks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

select.obj:    	select.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

voltra.obj:    	voltra.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bsstep.obj:    	bsstep.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fitexy.obj:    	fitexy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

laguer.obj:    	laguer.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

psdes.obj:     	psdes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

selip.obj:     	selip.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wt1.obj:       	wt1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

caldat.obj:    	caldat.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fixrts.obj:    	fixrts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lfit.obj:      	lfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pwt.obj:       	pwt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sfroid.obj:    	sfroid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wtn.obj:       	wtn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chder.obj:     	chder.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fleg.obj:      	fleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

linbcg.obj:    	linbcg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pwtset.obj:    	pwtset.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shell.obj:     	shell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wwghts.obj:    	wwghts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebev.obj:    	chebev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

flmoon.obj:    	flmoon.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

linmin.obj:    	linmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pythag.obj:    	pythag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shoot.obj:     	shoot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrac.obj:     	zbrac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebft.obj:    	chebft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fmin.obj:      	fmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lnsrch.obj:    	lnsrch.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pzextr.obj:    	pzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shootf.obj:    	shootf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrak.obj:     	zbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebpc.obj:    	chebpc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

four1.obj:     	four1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

locate.obj:    	locate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qgaus.obj:     	qgaus.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp1.obj:     	simp1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrent.obj:    	zbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chint.obj:     	chint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourew.obj:    	fourew.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lop.obj:       	lop.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrdcmp.obj:    	qrdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp2.obj:     	simp2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zrhqr.obj:     	zrhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chixy.obj:     	chixy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourfs.obj:    	fourfs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lubksb.obj:    	lubksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qromb.obj:     	qromb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp3.obj:     	simp3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zriddr.obj:    	zriddr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

choldc.obj:    	choldc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourn.obj:     	fourn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ludcmp.obj:    	ludcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qromo.obj:     	qromo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simplx.obj:    	simplx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zroots.obj:    	zroots.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cholsl.obj:    	cholsl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fpoly.obj:     	fpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

machar.obj:    	machar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qroot.obj:     	qroot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simpr.obj:     	simpr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chsone.obj:    	chsone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fred2.obj:     	fred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

matadd.obj:    	matadd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrsolv.obj:    	qrsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sinft.obj:     	sinft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chstwo.obj:    	chstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fredex.obj:    	fredex.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

matsub.obj:    	matsub.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrupdt.obj:    	qrupdt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

slvsm2.obj:    	slvsm2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cisi.obj:      	cisi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fredin.obj:    	fredin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

medfit.obj:    	medfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qsimp.obj:     	qsimp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

slvsml.obj:    	slvsml.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cntab1.obj:    	cntab1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

frenel.obj:    	frenel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

memcof.obj:    	memcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qtrap.obj:     	qtrap.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sncndn.obj:    	sncndn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsphfpt.obj:   	xsphfpt.cpp nrtypes.h nrtyp_nm.h newmat.h nrutil.h include.h myexcept.h nruti_nm.h

nruti_nm.obj:  	nruti_nm.cpp nr.h nruti_nm.h nrutil.h nrtypes.h nrtyp_nm.h newmat.h include.h myexcept.h

nrtyp_nm.obj:  	nrtyp_nm.cpp nr.h nrutil_nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

newmat1.obj:   	newmat1.cpp newmat.h include.h myexcept.h

newmat2.obj:   	newmat2.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat3.obj:   	newmat3.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat4.obj:   	newmat4.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat5.obj:   	newmat5.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat6.obj:   	newmat6.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat7.obj:   	newmat7.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat8.obj:   	newmat8.cpp include.h newmat.h newmatrc.h precisio.h myexcept.h controlw.h

newmatex.obj:  	newmatex.cpp include.h newmat.h myexcept.h

bandmat.obj:   	bandmat.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

submat.obj:    	submat.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

myexcept.obj:  	myexcept.cpp include.h myexcept.h

nrx_ex1.obj:   	nrx_ex1.cpp newmatap.h nr.h newmat.h include.h myexcept.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h

cholesky.obj:  	cholesky.cpp include.h newmat.h newmatrm.h myexcept.h

evalue.obj:    	evalue.cpp include.h newmatap.h newmatrm.h precisio.h newmat.h myexcept.h

fft.obj:       	fft.cpp include.h newmatap.h newmat.h myexcept.h

hholder.obj:   	hholder.cpp include.h newmatap.h newmat.h myexcept.h

jacobi.obj:    	jacobi.cpp include.h newmatap.h precisio.h newmatrm.h newmat.h myexcept.h

newfft.obj:    	newfft.cpp newmatap.h newmat.h include.h myexcept.h

sort.obj:      	sort.cpp include.h newmatap.h newmat.h myexcept.h

svd.obj:       	svd.cpp include.h newmatap.h newmatrm.h precisio.h newmat.h myexcept.h

newmatrm.obj:  	newmatrm.cpp newmat.h newmatrm.h include.h myexcept.h

xcyclic.obj:   	xcyclic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaulag.obj:   	xgaulag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpcshft.obj:   	xpcshft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xshell.obj:    	xshell.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xairy.obj:     	xairy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdawson.obj:   	xdawson.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgauleg.obj:   	xgauleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpearsn.obj:   	xpearsn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsimplx.obj:   	xsimplx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xamebsa.obj:   	xamebsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdbrent.obj:   	xdbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaussj.obj:   	xgaussj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xperiod.obj:   	xperiod.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsimpr.obj:    	xsimpr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xamoeba.obj:   	xamoeba.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xddpoly.obj:   	xddpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgcf.obj:      	xgcf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpiksr2.obj:   	xpiksr2.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsinft.obj:    	xsinft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xanneal.obj:   	xanneal.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdecchk.obj:   	xdecchk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgolden.obj:   	xgolden.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpiksrt.obj:   	xpiksrt.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsncndn.obj:   	xsncndn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xarcode.obj:   	xarcode.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdf1dim.obj:   	xdf1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgser.obj:     	xgser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xplgndr.obj:   	xplgndr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsobseq.obj:   	xsobseq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xavevar.obj:   	xavevar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdfpmin.obj:   	xdfpmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhpsel.obj:    	xhpsel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpoidev.obj:   	xpoidev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsor.obj:      	xsor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbalanc.obj:   	xbalanc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdfridr.obj:   	xdfridr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhpsort.obj:   	xhpsort.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolcoe.obj:   	xpolcoe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort.obj:     	xsort.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbandec.obj:   	xbandec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdftint.obj:   	xdftint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhqr.obj:      	xhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolcof.obj:   	xpolcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort2.obj:    	xsort2.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbanmul.obj:   	xbanmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdlinmin.obj:  	xdlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhuffman.obj:  	xhuffman.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpoldiv.obj:   	xpoldiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort3.obj:    	xsort3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbcucof.obj:   	xbcucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeclass.obj:   	xeclass.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhunt.obj:     	xhunt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolin2.obj:   	xpolin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspctrm.obj:   	xspctrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbcuint.obj:   	xbcuint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeclazz.obj:   	xeclazz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhypgeo.obj:   	xhypgeo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolint.obj:   	xpolint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspear.obj:    	xspear.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbeschb.obj:   	xbeschb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xei.obj:       	xei.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xicrc.obj:     	xicrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpowell.obj:   	xpowell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsphbes.obj:   	xsphbes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi.obj:    	xbessi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeigsrt.obj:   	xeigsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xigray.obj:    	xigray.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpredic.obj:   	xpredic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi0.obj:   	xbessi0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xelle.obj:     	xelle.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xindexx.obj:   	xindexx.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xprobks.obj:   	xprobks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplie2.obj:   	xsplie2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi1.obj:   	xbessi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xellf.obj:     	xellf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xirbit1.obj:   	xirbit1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpsdes.obj:    	xpsdes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplin2.obj:   	xsplin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessik.obj:   	xbessik.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xellpi.obj:    	xellpi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xirbit2.obj:   	xirbit2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpzextr.obj:   	xpzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspline.obj:   	xspline.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj.obj:    	xbessj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xelmhes.obj:   	xelmhes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xjacobi.obj:   	xjacobi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqgaus.obj:    	xqgaus.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplint.obj:   	xsplint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj0.obj:   	xbessj0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerfcc.obj:    	xerfcc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xjulday.obj:   	xjulday.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrdcmp.obj:   	xqrdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprsax.obj:   	xsprsax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj1.obj:   	xbessj1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerff.obj:     	xerff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkendl1.obj:   	xkendl1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqromb.obj:    	xqromb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprsin.obj:   	xsprsin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessjy.obj:   	xbessjy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerffc.obj:    	xerffc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkendl2.obj:   	xkendl2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqromo.obj:    	xqromo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprspm.obj:   	xsprspm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk.obj:    	xbessk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeulsum.obj:   	xeulsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xks2d1s.obj:   	xks2d1s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqroot.obj:    	xqroot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstm.obj:   	xsprstm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk0.obj:   	xbessk0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xevlmem.obj:   	xevlmem.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xks2d2s.obj:   	xks2d2s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrsolv.obj:   	xqrsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstp.obj:   	xsprstp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk1.obj:   	xbessk1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xexpdev.obj:   	xexpdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xksone.obj:    	xksone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrupdt.obj:   	xqrupdt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstx.obj:   	xsprstx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy.obj:    	xbessy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xexpint.obj:   	xexpint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkstwo.obj:    	xkstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqsimp.obj:    	xqsimp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstifbs.obj:   	xstifbs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy0.obj:   	xbessy0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xf1dim.obj:    	xf1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlaguer.obj:   	xlaguer.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqtrap.obj:    	xqtrap.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstiff.obj:    	xstiff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy1.obj:   	xbessy1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfactln.obj:   	xfactln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlfit.obj:     	xlfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xquad3d.obj:   	xquad3d.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstoerm.obj:   	xstoerm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbeta.obj:     	xbeta.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfactrl.obj:   	xfactrl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlinbcg.obj:   	xlinbcg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xran.obj:      	xran.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvbksb.obj:   	xsvbksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbetai.obj:    	xbetai.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfasper.obj:   	xfasper.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlinmin.obj:   	xlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xran4.obj:     	xran4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdcmp.obj:   	xsvdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbico.obj:     	xbico.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfgauss.obj:   	xfgauss.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlocate.obj:   	xlocate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrank.obj:     	xrank.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdfit.obj:   	xsvdfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbnldev.obj:   	xbnldev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfit.obj:      	xfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlubksb.obj:   	xlubksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xratint.obj:   	xratint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdvar.obj:   	xsvdvar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbrent.obj:    	xbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfitexy.obj:   	xfitexy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xludcmp.obj:   	xludcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xratlsq.obj:   	xratlsq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtoeplz.obj:   	xtoeplz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbroydn.obj:   	xbroydn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfixrts.obj:   	xfixrts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmachar.obj:   	xmachar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrc.obj:       	xrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtptest.obj:   	xtptest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbsstep.obj:   	xbsstep.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfleg.obj:     	xfleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmedfit.obj:   	xmedfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrd.obj:       	xrd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtqli.obj:     	xtqli.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcaldat.obj:   	xcaldat.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xflmoon.obj:   	xflmoon.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmemcof.obj:   	xmemcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrealft.obj:   	xrealft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtrapzd.obj:   	xtrapzd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchder.obj:    	xchder.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfour1.obj:    	xfour1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmgfas.obj:    	xmgfas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrf.obj:       	xrf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtred2.obj:    	xtred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebev.obj:   	xchebev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfourfs.obj:   	xfourfs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmglin.obj:    	xmglin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrj.obj:       	xrj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtridag.obj:   	xtridag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebft.obj:   	xchebft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfourn.obj:    	xfourn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmidpnt.obj:   	xmidpnt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrk4.obj:      	xrk4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xttest.obj:    	xttest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebpc.obj:   	xchebpc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfpoly.obj:    	xfpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmiser.obj:    	xmiser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrkdumb.obj:   	xrkdumb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtutest.obj:   	xtutest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchint.obj:    	xchint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfred2.obj:    	xfred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmmid.obj:     	xmmid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrkqs.obj:     	xrkqs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtwofft.obj:   	xtwofft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcholsl.obj:   	xcholsl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfredin.obj:   	xfredin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmnbrak.obj:   	xmnbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrlft3.obj:    	xrlft3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvander.obj:   	xvander.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchsone.obj:   	xchsone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfrenel.obj:   	xfrenel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmnewt.obj:    	xmnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrofunc.obj:   	xrofunc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvegas.obj:    	xvegas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchstwo.obj:   	xchstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfrprmn.obj:   	xfrprmn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmoment.obj:   	xmoment.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtbis.obj:    	xrtbis.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvoltra.obj:   	xvoltra.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcisi.obj:     	xcisi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xftest.obj:    	xftest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmppi.obj:     	xmppi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtflsp.obj:   	xrtflsp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xwt1.obj:      	xwt1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcntab1.obj:   	xcntab1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgamdev.obj:   	xgamdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmprove.obj:   	xmprove.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtnewt.obj:   	xrtnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xwtn.obj:      	xwtn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcntab2.obj:   	xcntab2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammln.obj:   	xgammln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmrqcof.obj:   	xmrqcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtsafe.obj:   	xrtsafe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrac.obj:    	xzbrac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xconvlv.obj:   	xconvlv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammp.obj:    	xgammp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmrqmin.obj:   	xmrqmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtsec.obj:    	xrtsec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrak.obj:    	xzbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcorrel.obj:   	xcorrel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammq.obj:    	xgammq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xnewt.obj:     	xnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrzextr.obj:   	xrzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrent.obj:   	xzbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcosft1.obj:   	xcosft1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgasdev.obj:   	xgasdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xodeint.obj:   	xodeint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsavgol.obj:   	xsavgol.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzrhqr.obj:    	xzrhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcosft2.obj:   	xcosft2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaucof.obj:   	xgaucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xorthog.obj:   	xorthog.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xscrsho.obj:   	xscrsho.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzriddr.obj:   	xzriddr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcovsrt.obj:   	xcovsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgauher.obj:   	xgauher.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpade.obj:     	xpade.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xselect.obj:   	xselect.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzroots.obj:   	xzroots.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcrank.obj:    	xcrank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaujac.obj:   	xgaujac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpccheb.obj:   	xpccheb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xselip.obj:    	xselip.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

nrx_ex1.txx:   	nrx_ex1.exe
		$(PRE)nrx_ex1 > nrx_ex1.txx
		$(DIFF) nrx_ex1.txt nrx_ex1.txx

xcyclic.txx:   	xcyclic.exe
		$(PRE)xcyclic > xcyclic.txx
		$(DIFF) xcyclic.txt xcyclic.txx

xgaulag.txx:   	xgaulag.exe
		$(PRE)xgaulag > xgaulag.txx
		$(DIFF) xgaulag.txt xgaulag.txx

xpcshft.txx:   	xpcshft.exe
		$(PRE)xpcshft > xpcshft.txx
		$(DIFF) xpcshft.txt xpcshft.txx

xshell.txx:    	xshell.exe
		$(PRE)xshell > xshell.txx
		$(DIFF) xshell.txt xshell.txx

xairy.txx:     	xairy.exe
		$(PRE)xairy > xairy.txx
		$(DIFF) xairy.txt xairy.txx

xdawson.txx:   	xdawson.exe
		$(PRE)xdawson > xdawson.txx
		$(DIFF) xdawson.txt xdawson.txx

xgauleg.txx:   	xgauleg.exe
		$(PRE)xgauleg > xgauleg.txx
		$(DIFF) xgauleg.txt xgauleg.txx

xpearsn.txx:   	xpearsn.exe
		$(PRE)xpearsn > xpearsn.txx
		$(DIFF) xpearsn.txt xpearsn.txx

xsimplx.txx:   	xsimplx.exe
		$(PRE)xsimplx > xsimplx.txx
		$(DIFF) xsimplx.txt xsimplx.txx

xamebsa.txx:   	xamebsa.exe
		$(PRE)xamebsa > xamebsa.txx
		$(DIFF) xamebsa.txt xamebsa.txx

xdbrent.txx:   	xdbrent.exe
		$(PRE)xdbrent > xdbrent.txx
		$(DIFF) xdbrent.txt xdbrent.txx

xgaussj.txx:   	xgaussj.exe
		$(PRE)xgaussj > xgaussj.txx
		$(DIFF) xgaussj.txt xgaussj.txx

xperiod.txx:   	xperiod.exe
		$(PRE)xperiod > xperiod.txx
		$(DIFF) xperiod.txt xperiod.txx

xsimpr.txx:    	xsimpr.exe
		$(PRE)xsimpr > xsimpr.txx
		$(DIFF) xsimpr.txt xsimpr.txx

xamoeba.txx:   	xamoeba.exe
		$(PRE)xamoeba > xamoeba.txx
		$(DIFF) xamoeba.txt xamoeba.txx

xddpoly.txx:   	xddpoly.exe
		$(PRE)xddpoly > xddpoly.txx
		$(DIFF) xddpoly.txt xddpoly.txx

xgcf.txx:      	xgcf.exe
		$(PRE)xgcf > xgcf.txx
		$(DIFF) xgcf.txt xgcf.txx

xpiksr2.txx:   	xpiksr2.exe
		$(PRE)xpiksr2 > xpiksr2.txx
		$(DIFF) xpiksr2.txt xpiksr2.txx

xsinft.txx:    	xsinft.exe
		$(PRE)xsinft > xsinft.txx
		$(DIFF) xsinft.txt xsinft.txx

xanneal.txx:   	xanneal.exe
		$(PRE)xanneal > xanneal.txx
		$(DIFF) xanneal.txt xanneal.txx

xdecchk.txx:   	xdecchk.exe
		$(PRE)xdecchk > xdecchk.txx
		$(DIFF) xdecchk.txt xdecchk.txx

xgolden.txx:   	xgolden.exe
		$(PRE)xgolden > xgolden.txx
		$(DIFF) xgolden.txt xgolden.txx

xpiksrt.txx:   	xpiksrt.exe
		$(PRE)xpiksrt > xpiksrt.txx
		$(DIFF) xpiksrt.txt xpiksrt.txx

xsncndn.txx:   	xsncndn.exe
		$(PRE)xsncndn > xsncndn.txx
		$(DIFF) xsncndn.txt xsncndn.txx

xarcode.txx:   	xarcode.exe
		$(PRE)xarcode > xarcode.txx
		$(DIFF) xarcode.txt xarcode.txx

xdf1dim.txx:   	xdf1dim.exe
		$(PRE)xdf1dim > xdf1dim.txx
		$(DIFF) xdf1dim.txt xdf1dim.txx

xgser.txx:     	xgser.exe
		$(PRE)xgser > xgser.txx
		$(DIFF) xgser.txt xgser.txx

xplgndr.txx:   	xplgndr.exe
		$(PRE)xplgndr > xplgndr.txx
		$(DIFF) xplgndr.txt xplgndr.txx

xsobseq.txx:   	xsobseq.exe
		$(PRE)xsobseq > xsobseq.txx
		$(DIFF) xsobseq.txt xsobseq.txx

xavevar.txx:   	xavevar.exe
		$(PRE)xavevar > xavevar.txx
		$(DIFF) xavevar.txt xavevar.txx

xdfpmin.txx:   	xdfpmin.exe
		$(PRE)xdfpmin > xdfpmin.txx
		$(DIFF) xdfpmin.txt xdfpmin.txx

xhpsel.txx:    	xhpsel.exe
		$(PRE)xhpsel > xhpsel.txx
		$(DIFF) xhpsel.txt xhpsel.txx

xpoidev.txx:   	xpoidev.exe
		$(PRE)xpoidev > xpoidev.txx
		$(DIFF) xpoidev.txt xpoidev.txx

xsor.txx:      	xsor.exe
		$(PRE)xsor > xsor.txx
		$(DIFF) xsor.txt xsor.txx

xbalanc.txx:   	xbalanc.exe
		$(PRE)xbalanc > xbalanc.txx
		$(DIFF) xbalanc.txt xbalanc.txx

xdfridr.txx:   	xdfridr.exe
		$(PRE)xdfridr > xdfridr.txx
		$(DIFF) xdfridr.txt xdfridr.txx

xhpsort.txx:   	xhpsort.exe
		$(PRE)xhpsort > xhpsort.txx
		$(DIFF) xhpsort.txt xhpsort.txx

xpolcoe.txx:   	xpolcoe.exe
		$(PRE)xpolcoe > xpolcoe.txx
		$(DIFF) xpolcoe.txt xpolcoe.txx

xsort.txx:     	xsort.exe
		$(PRE)xsort > xsort.txx
		$(DIFF) xsort.txt xsort.txx

xbandec.txx:   	xbandec.exe
		$(PRE)xbandec > xbandec.txx
		$(DIFF) xbandec.txt xbandec.txx

xdftint.txx:   	xdftint.exe
		$(PRE)xdftint > xdftint.txx
		$(DIFF) xdftint.txt xdftint.txx

xhqr.txx:      	xhqr.exe
		$(PRE)xhqr > xhqr.txx
		$(DIFF) xhqr.txt xhqr.txx

xpolcof.txx:   	xpolcof.exe
		$(PRE)xpolcof > xpolcof.txx
		$(DIFF) xpolcof.txt xpolcof.txx

xsort2.txx:    	xsort2.exe
		$(PRE)xsort2 > xsort2.txx
		$(DIFF) xsort2.txt xsort2.txx

xbanmul.txx:   	xbanmul.exe
		$(PRE)xbanmul > xbanmul.txx
		$(DIFF) xbanmul.txt xbanmul.txx

xdlinmin.txx:  	xdlinmin.exe
		$(PRE)xdlinmin > xdlinmin.txx
		$(DIFF) xdlinmin.txt xdlinmin.txx

xhuffman.txx:  	xhuffman.exe
		$(PRE)xhuffman > xhuffman.txx
		$(DIFF) xhuffman.txt xhuffman.txx

xpoldiv.txx:   	xpoldiv.exe
		$(PRE)xpoldiv > xpoldiv.txx
		$(DIFF) xpoldiv.txt xpoldiv.txx

xsort3.txx:    	xsort3.exe
		$(PRE)xsort3 > xsort3.txx
		$(DIFF) xsort3.txt xsort3.txx

xbcucof.txx:   	xbcucof.exe
		$(PRE)xbcucof > xbcucof.txx
		$(DIFF) xbcucof.txt xbcucof.txx

xeclass.txx:   	xeclass.exe
		$(PRE)xeclass > xeclass.txx
		$(DIFF) xeclass.txt xeclass.txx

xhunt.txx:     	xhunt.exe
		$(PRE)xhunt > xhunt.txx
		$(DIFF) xhunt.txt xhunt.txx

xpolin2.txx:   	xpolin2.exe
		$(PRE)xpolin2 > xpolin2.txx
		$(DIFF) xpolin2.txt xpolin2.txx

xspctrm.txx:   	xspctrm.exe
		$(PRE)xspctrm > xspctrm.txx
		$(DIFF) xspctrm.txt xspctrm.txx

xbcuint.txx:   	xbcuint.exe
		$(PRE)xbcuint > xbcuint.txx
		$(DIFF) xbcuint.txt xbcuint.txx

xeclazz.txx:   	xeclazz.exe
		$(PRE)xeclazz > xeclazz.txx
		$(DIFF) xeclazz.txt xeclazz.txx

xhypgeo.txx:   	xhypgeo.exe
		$(PRE)xhypgeo > xhypgeo.txx
		$(DIFF) xhypgeo.txt xhypgeo.txx

xpolint.txx:   	xpolint.exe
		$(PRE)xpolint > xpolint.txx
		$(DIFF) xpolint.txt xpolint.txx

xspear.txx:    	xspear.exe
		$(PRE)xspear > xspear.txx
		$(DIFF) xspear.txt xspear.txx

xbeschb.txx:   	xbeschb.exe
		$(PRE)xbeschb > xbeschb.txx
		$(DIFF) xbeschb.txt xbeschb.txx

xei.txx:       	xei.exe
		$(PRE)xei > xei.txx
		$(DIFF) xei.txt xei.txx

xicrc.txx:     	xicrc.exe
		$(PRE)xicrc > xicrc.txx
		$(DIFF) xicrc.txt xicrc.txx

xpowell.txx:   	xpowell.exe
		$(PRE)xpowell > xpowell.txx
		$(DIFF) xpowell.txt xpowell.txx

xsphbes.txx:   	xsphbes.exe
		$(PRE)xsphbes > xsphbes.txx
		$(DIFF) xsphbes.txt xsphbes.txx

xbessi.txx:    	xbessi.exe
		$(PRE)xbessi > xbessi.txx
		$(DIFF) xbessi.txt xbessi.txx

xeigsrt.txx:   	xeigsrt.exe
		$(PRE)xeigsrt > xeigsrt.txx
		$(DIFF) xeigsrt.txt xeigsrt.txx

xigray.txx:    	xigray.exe
		$(PRE)xigray > xigray.txx
		$(DIFF) xigray.txt xigray.txx

xpredic.txx:   	xpredic.exe
		$(PRE)xpredic > xpredic.txx
		$(DIFF) xpredic.txt xpredic.txx

xbessi0.txx:   	xbessi0.exe
		$(PRE)xbessi0 > xbessi0.txx
		$(DIFF) xbessi0.txt xbessi0.txx

xelle.txx:     	xelle.exe
		$(PRE)xelle > xelle.txx
		$(DIFF) xelle.txt xelle.txx

xindexx.txx:   	xindexx.exe
		$(PRE)xindexx > xindexx.txx
		$(DIFF) xindexx.txt xindexx.txx

xprobks.txx:   	xprobks.exe
		$(PRE)xprobks > xprobks.txx
		$(DIFF) xprobks.txt xprobks.txx

xsplie2.txx:   	xsplie2.exe
		$(PRE)xsplie2 > xsplie2.txx
		$(DIFF) xsplie2.txt xsplie2.txx

xbessi1.txx:   	xbessi1.exe
		$(PRE)xbessi1 > xbessi1.txx
		$(DIFF) xbessi1.txt xbessi1.txx

xellf.txx:     	xellf.exe
		$(PRE)xellf > xellf.txx
		$(DIFF) xellf.txt xellf.txx

xirbit1.txx:   	xirbit1.exe
		$(PRE)xirbit1 > xirbit1.txx
		$(DIFF) xirbit1.txt xirbit1.txx

xpsdes.txx:    	xpsdes.exe
		$(PRE)xpsdes > xpsdes.txx
		$(DIFF) xpsdes.txt xpsdes.txx

xsplin2.txx:   	xsplin2.exe
		$(PRE)xsplin2 > xsplin2.txx
		$(DIFF) xsplin2.txt xsplin2.txx

xbessik.txx:   	xbessik.exe
		$(PRE)xbessik > xbessik.txx
		$(DIFF) xbessik.txt xbessik.txx

xellpi.txx:    	xellpi.exe
		$(PRE)xellpi > xellpi.txx
		$(DIFF) xellpi.txt xellpi.txx

xirbit2.txx:   	xirbit2.exe
		$(PRE)xirbit2 > xirbit2.txx
		$(DIFF) xirbit2.txt xirbit2.txx

xpzextr.txx:   	xpzextr.exe
		$(PRE)xpzextr > xpzextr.txx
		$(DIFF) xpzextr.txt xpzextr.txx

xspline.txx:   	xspline.exe
		$(PRE)xspline > xspline.txx
		$(DIFF) xspline.txt xspline.txx

xbessj.txx:    	xbessj.exe
		$(PRE)xbessj > xbessj.txx
		$(DIFF) xbessj.txt xbessj.txx

xelmhes.txx:   	xelmhes.exe
		$(PRE)xelmhes > xelmhes.txx
		$(DIFF) xelmhes.txt xelmhes.txx

xjacobi.txx:   	xjacobi.exe
		$(PRE)xjacobi > xjacobi.txx
		$(DIFF) xjacobi.txt xjacobi.txx

xqgaus.txx:    	xqgaus.exe
		$(PRE)xqgaus > xqgaus.txx
		$(DIFF) xqgaus.txt xqgaus.txx

xsplint.txx:   	xsplint.exe
		$(PRE)xsplint > xsplint.txx
		$(DIFF) xsplint.txt xsplint.txx

xbessj0.txx:   	xbessj0.exe
		$(PRE)xbessj0 > xbessj0.txx
		$(DIFF) xbessj0.txt xbessj0.txx

xerfcc.txx:    	xerfcc.exe
		$(PRE)xerfcc > xerfcc.txx
		$(DIFF) xerfcc.txt xerfcc.txx

xjulday.txx:   	xjulday.exe
		$(PRE)xjulday > xjulday.txx
		$(DIFF) xjulday.txt xjulday.txx

xqrdcmp.txx:   	xqrdcmp.exe
		$(PRE)xqrdcmp > xqrdcmp.txx
		$(DIFF) xqrdcmp.txt xqrdcmp.txx

xsprsax.txx:   	xsprsax.exe
		$(PRE)xsprsax > xsprsax.txx
		$(DIFF) xsprsax.txt xsprsax.txx

xbessj1.txx:   	xbessj1.exe
		$(PRE)xbessj1 > xbessj1.txx
		$(DIFF) xbessj1.txt xbessj1.txx

xerff.txx:     	xerff.exe
		$(PRE)xerff > xerff.txx
		$(DIFF) xerff.txt xerff.txx

xkendl1.txx:   	xkendl1.exe
		$(PRE)xkendl1 > xkendl1.txx
		$(DIFF) xkendl1.txt xkendl1.txx

xqromb.txx:    	xqromb.exe
		$(PRE)xqromb > xqromb.txx
		$(DIFF) xqromb.txt xqromb.txx

xsprsin.txx:   	xsprsin.exe
		$(PRE)xsprsin > xsprsin.txx
		$(DIFF) xsprsin.txt xsprsin.txx

xbessjy.txx:   	xbessjy.exe
		$(PRE)xbessjy > xbessjy.txx
		$(DIFF) xbessjy.txt xbessjy.txx

xerffc.txx:    	xerffc.exe
		$(PRE)xerffc > xerffc.txx
		$(DIFF) xerffc.txt xerffc.txx

xkendl2.txx:   	xkendl2.exe
		$(PRE)xkendl2 > xkendl2.txx
		$(DIFF) xkendl2.txt xkendl2.txx

xqromo.txx:    	xqromo.exe
		$(PRE)xqromo > xqromo.txx
		$(DIFF) xqromo.txt xqromo.txx

xsprspm.txx:   	xsprspm.exe
		$(PRE)xsprspm > xsprspm.txx
		$(DIFF) xsprspm.txt xsprspm.txx

xbessk.txx:    	xbessk.exe
		$(PRE)xbessk > xbessk.txx
		$(DIFF) xbessk.txt xbessk.txx

xeulsum.txx:   	xeulsum.exe
		$(PRE)xeulsum > xeulsum.txx
		$(DIFF) xeulsum.txt xeulsum.txx

xks2d1s.txx:   	xks2d1s.exe
		$(PRE)xks2d1s > xks2d1s.txx
		$(DIFF) xks2d1s.txt xks2d1s.txx

xqroot.txx:    	xqroot.exe
		$(PRE)xqroot > xqroot.txx
		$(DIFF) xqroot.txt xqroot.txx

xsprstm.txx:   	xsprstm.exe
		$(PRE)xsprstm > xsprstm.txx
		$(DIFF) xsprstm.txt xsprstm.txx

xbessk0.txx:   	xbessk0.exe
		$(PRE)xbessk0 > xbessk0.txx
		$(DIFF) xbessk0.txt xbessk0.txx

xevlmem.txx:   	xevlmem.exe
		$(PRE)xevlmem > xevlmem.txx
		$(DIFF) xevlmem.txt xevlmem.txx

xks2d2s.txx:   	xks2d2s.exe
		$(PRE)xks2d2s > xks2d2s.txx
		$(DIFF) xks2d2s.txt xks2d2s.txx

xqrsolv.txx:   	xqrsolv.exe
		$(PRE)xqrsolv > xqrsolv.txx
		$(DIFF) xqrsolv.txt xqrsolv.txx

xsprstp.txx:   	xsprstp.exe
		$(PRE)xsprstp > xsprstp.txx
		$(DIFF) xsprstp.txt xsprstp.txx

xbessk1.txx:   	xbessk1.exe
		$(PRE)xbessk1 > xbessk1.txx
		$(DIFF) xbessk1.txt xbessk1.txx

xexpdev.txx:   	xexpdev.exe
		$(PRE)xexpdev > xexpdev.txx
		$(DIFF) xexpdev.txt xexpdev.txx

xksone.txx:    	xksone.exe
		$(PRE)xksone > xksone.txx
		$(DIFF) xksone.txt xksone.txx

xqrupdt.txx:   	xqrupdt.exe
		$(PRE)xqrupdt > xqrupdt.txx
		$(DIFF) xqrupdt.txt xqrupdt.txx

xsprstx.txx:   	xsprstx.exe
		$(PRE)xsprstx > xsprstx.txx
		$(DIFF) xsprstx.txt xsprstx.txx

xbessy.txx:    	xbessy.exe
		$(PRE)xbessy > xbessy.txx
		$(DIFF) xbessy.txt xbessy.txx

xexpint.txx:   	xexpint.exe
		$(PRE)xexpint > xexpint.txx
		$(DIFF) xexpint.txt xexpint.txx

xkstwo.txx:    	xkstwo.exe
		$(PRE)xkstwo > xkstwo.txx
		$(DIFF) xkstwo.txt xkstwo.txx

xqsimp.txx:    	xqsimp.exe
		$(PRE)xqsimp > xqsimp.txx
		$(DIFF) xqsimp.txt xqsimp.txx

xstifbs.txx:   	xstifbs.exe
		$(PRE)xstifbs > xstifbs.txx
		$(DIFF) xstifbs.txt xstifbs.txx

xbessy0.txx:   	xbessy0.exe
		$(PRE)xbessy0 > xbessy0.txx
		$(DIFF) xbessy0.txt xbessy0.txx

xf1dim.txx:    	xf1dim.exe
		$(PRE)xf1dim > xf1dim.txx
		$(DIFF) xf1dim.txt xf1dim.txx

xlaguer.txx:   	xlaguer.exe
		$(PRE)xlaguer > xlaguer.txx
		$(DIFF) xlaguer.txt xlaguer.txx

xqtrap.txx:    	xqtrap.exe
		$(PRE)xqtrap > xqtrap.txx
		$(DIFF) xqtrap.txt xqtrap.txx

xstiff.txx:    	xstiff.exe
		$(PRE)xstiff > xstiff.txx
		$(DIFF) xstiff.txt xstiff.txx

xbessy1.txx:   	xbessy1.exe
		$(PRE)xbessy1 > xbessy1.txx
		$(DIFF) xbessy1.txt xbessy1.txx

xfactln.txx:   	xfactln.exe
		$(PRE)xfactln > xfactln.txx
		$(DIFF) xfactln.txt xfactln.txx

xlfit.txx:     	xlfit.exe
		$(PRE)xlfit > xlfit.txx
		$(DIFF) xlfit.txt xlfit.txx

xquad3d.txx:   	xquad3d.exe
		$(PRE)xquad3d > xquad3d.txx
		$(DIFF) xquad3d.txt xquad3d.txx

xstoerm.txx:   	xstoerm.exe
		$(PRE)xstoerm > xstoerm.txx
		$(DIFF) xstoerm.txt xstoerm.txx

xbeta.txx:     	xbeta.exe
		$(PRE)xbeta > xbeta.txx
		$(DIFF) xbeta.txt xbeta.txx

xfactrl.txx:   	xfactrl.exe
		$(PRE)xfactrl > xfactrl.txx
		$(DIFF) xfactrl.txt xfactrl.txx

xlinbcg.txx:   	xlinbcg.exe
		$(PRE)xlinbcg > xlinbcg.txx
		$(DIFF) xlinbcg.txt xlinbcg.txx

xran.txx:      	xran.exe
		$(PRE)xran > xran.txx
		$(DIFF) xran.txt xran.txx

xsvbksb.txx:   	xsvbksb.exe
		$(PRE)xsvbksb > xsvbksb.txx
		$(DIFF) xsvbksb.txt xsvbksb.txx

xbetai.txx:    	xbetai.exe
		$(PRE)xbetai > xbetai.txx
		$(DIFF) xbetai.txt xbetai.txx

xfasper.txx:   	xfasper.exe
		$(PRE)xfasper > xfasper.txx
		$(DIFF) xfasper.txt xfasper.txx

xlinmin.txx:   	xlinmin.exe
		$(PRE)xlinmin > xlinmin.txx
		$(DIFF) xlinmin.txt xlinmin.txx

xran4.txx:     	xran4.exe
		$(PRE)xran4 > xran4.txx
		$(DIFF) xran4.txt xran4.txx

xsvdcmp.txx:   	xsvdcmp.exe
		$(PRE)xsvdcmp > xsvdcmp.txx
		$(DIFF) xsvdcmp.txt xsvdcmp.txx

xbico.txx:     	xbico.exe
		$(PRE)xbico > xbico.txx
		$(DIFF) xbico.txt xbico.txx

xfgauss.txx:   	xfgauss.exe
		$(PRE)xfgauss > xfgauss.txx
		$(DIFF) xfgauss.txt xfgauss.txx

xlocate.txx:   	xlocate.exe
		$(PRE)xlocate > xlocate.txx
		$(DIFF) xlocate.txt xlocate.txx

xrank.txx:     	xrank.exe
		$(PRE)xrank > xrank.txx
		$(DIFF) xrank.txt xrank.txx

xsvdfit.txx:   	xsvdfit.exe
		$(PRE)xsvdfit > xsvdfit.txx
		$(DIFF) xsvdfit.txt xsvdfit.txx

xbnldev.txx:   	xbnldev.exe
		$(PRE)xbnldev > xbnldev.txx
		$(DIFF) xbnldev.txt xbnldev.txx

xfit.txx:      	xfit.exe
		$(PRE)xfit > xfit.txx
		$(DIFF) xfit.txt xfit.txx

xlubksb.txx:   	xlubksb.exe
		$(PRE)xlubksb > xlubksb.txx
		$(DIFF) xlubksb.txt xlubksb.txx

xratint.txx:   	xratint.exe
		$(PRE)xratint > xratint.txx
		$(DIFF) xratint.txt xratint.txx

xsvdvar.txx:   	xsvdvar.exe
		$(PRE)xsvdvar > xsvdvar.txx
		$(DIFF) xsvdvar.txt xsvdvar.txx

xbrent.txx:    	xbrent.exe
		$(PRE)xbrent > xbrent.txx
		$(DIFF) xbrent.txt xbrent.txx

xfitexy.txx:   	xfitexy.exe
		$(PRE)xfitexy > xfitexy.txx
		$(DIFF) xfitexy.txt xfitexy.txx

xludcmp.txx:   	xludcmp.exe
		$(PRE)xludcmp > xludcmp.txx
		$(DIFF) xludcmp.txt xludcmp.txx

xratlsq.txx:   	xratlsq.exe
		$(PRE)xratlsq > xratlsq.txx
		$(DIFF) xratlsq.txt xratlsq.txx

xtoeplz.txx:   	xtoeplz.exe
		$(PRE)xtoeplz > xtoeplz.txx
		$(DIFF) xtoeplz.txt xtoeplz.txx

xbroydn.txx:   	xbroydn.exe
		$(PRE)xbroydn > xbroydn.txx
		$(DIFF) xbroydn.txt xbroydn.txx

xfixrts.txx:   	xfixrts.exe
		$(PRE)xfixrts > xfixrts.txx
		$(DIFF) xfixrts.txt xfixrts.txx

xmachar.txx:   	xmachar.exe
		$(PRE)xmachar > xmachar.txx
		$(DIFF) xmachar.txt xmachar.txx

xrc.txx:       	xrc.exe
		$(PRE)xrc > xrc.txx
		$(DIFF) xrc.txt xrc.txx

xtptest.txx:   	xtptest.exe
		$(PRE)xtptest > xtptest.txx
		$(DIFF) xtptest.txt xtptest.txx

xbsstep.txx:   	xbsstep.exe
		$(PRE)xbsstep > xbsstep.txx
		$(DIFF) xbsstep.txt xbsstep.txx

xfleg.txx:     	xfleg.exe
		$(PRE)xfleg > xfleg.txx
		$(DIFF) xfleg.txt xfleg.txx

xmedfit.txx:   	xmedfit.exe
		$(PRE)xmedfit > xmedfit.txx
		$(DIFF) xmedfit.txt xmedfit.txx

xrd.txx:       	xrd.exe
		$(PRE)xrd > xrd.txx
		$(DIFF) xrd.txt xrd.txx

xtqli.txx:     	xtqli.exe
		$(PRE)xtqli > xtqli.txx
		$(DIFF) xtqli.txt xtqli.txx

xcaldat.txx:   	xcaldat.exe
		$(PRE)xcaldat > xcaldat.txx
		$(DIFF) xcaldat.txt xcaldat.txx

xflmoon.txx:   	xflmoon.exe
		$(PRE)xflmoon > xflmoon.txx
		$(DIFF) xflmoon.txt xflmoon.txx

xmemcof.txx:   	xmemcof.exe
		$(PRE)xmemcof > xmemcof.txx
		$(DIFF) xmemcof.txt xmemcof.txx

xrealft.txx:   	xrealft.exe
		$(PRE)xrealft > xrealft.txx
		$(DIFF) xrealft.txt xrealft.txx

xtrapzd.txx:   	xtrapzd.exe
		$(PRE)xtrapzd > xtrapzd.txx
		$(DIFF) xtrapzd.txt xtrapzd.txx

xchder.txx:    	xchder.exe
		$(PRE)xchder > xchder.txx
		$(DIFF) xchder.txt xchder.txx

xfour1.txx:    	xfour1.exe
		$(PRE)xfour1 > xfour1.txx
		$(DIFF) xfour1.txt xfour1.txx

xmgfas.txx:    	xmgfas.exe
		$(PRE)xmgfas > xmgfas.txx
		$(DIFF) xmgfas.txt xmgfas.txx

xrf.txx:       	xrf.exe
		$(PRE)xrf > xrf.txx
		$(DIFF) xrf.txt xrf.txx

xtred2.txx:    	xtred2.exe
		$(PRE)xtred2 > xtred2.txx
		$(DIFF) xtred2.txt xtred2.txx

xchebev.txx:   	xchebev.exe
		$(PRE)xchebev > xchebev.txx
		$(DIFF) xchebev.txt xchebev.txx

xfourfs.txx:   	xfourfs.exe
		$(PRE)xfourfs > xfourfs.txx
		$(DIFF) xfourfs.txt xfourfs.txx

xmglin.txx:    	xmglin.exe
		$(PRE)xmglin > xmglin.txx
		$(DIFF) xmglin.txt xmglin.txx

xrj.txx:       	xrj.exe
		$(PRE)xrj > xrj.txx
		$(DIFF) xrj.txt xrj.txx

xtridag.txx:   	xtridag.exe
		$(PRE)xtridag > xtridag.txx
		$(DIFF) xtridag.txt xtridag.txx

xchebft.txx:   	xchebft.exe
		$(PRE)xchebft > xchebft.txx
		$(DIFF) xchebft.txt xchebft.txx

xfourn.txx:    	xfourn.exe
		$(PRE)xfourn > xfourn.txx
		$(DIFF) xfourn.txt xfourn.txx

xmidpnt.txx:   	xmidpnt.exe
		$(PRE)xmidpnt > xmidpnt.txx
		$(DIFF) xmidpnt.txt xmidpnt.txx

xrk4.txx:      	xrk4.exe
		$(PRE)xrk4 > xrk4.txx
		$(DIFF) xrk4.txt xrk4.txx

xttest.txx:    	xttest.exe
		$(PRE)xttest > xttest.txx
		$(DIFF) xttest.txt xttest.txx

xchebpc.txx:   	xchebpc.exe
		$(PRE)xchebpc > xchebpc.txx
		$(DIFF) xchebpc.txt xchebpc.txx

xfpoly.txx:    	xfpoly.exe
		$(PRE)xfpoly > xfpoly.txx
		$(DIFF) xfpoly.txt xfpoly.txx

xmiser.txx:    	xmiser.exe
		$(PRE)xmiser > xmiser.txx
		$(DIFF) xmiser.txt xmiser.txx

xrkdumb.txx:   	xrkdumb.exe
		$(PRE)xrkdumb > xrkdumb.txx
		$(DIFF) xrkdumb.txt xrkdumb.txx

xtutest.txx:   	xtutest.exe
		$(PRE)xtutest > xtutest.txx
		$(DIFF) xtutest.txt xtutest.txx

xchint.txx:    	xchint.exe
		$(PRE)xchint > xchint.txx
		$(DIFF) xchint.txt xchint.txx

xfred2.txx:    	xfred2.exe
		$(PRE)xfred2 > xfred2.txx
		$(DIFF) xfred2.txt xfred2.txx

xmmid.txx:     	xmmid.exe
		$(PRE)xmmid > xmmid.txx
		$(DIFF) xmmid.txt xmmid.txx

xrkqs.txx:     	xrkqs.exe
		$(PRE)xrkqs > xrkqs.txx
		$(DIFF) xrkqs.txt xrkqs.txx

xtwofft.txx:   	xtwofft.exe
		$(PRE)xtwofft > xtwofft.txx
		$(DIFF) xtwofft.txt xtwofft.txx

xcholsl.txx:   	xcholsl.exe
		$(PRE)xcholsl > xcholsl.txx
		$(DIFF) xcholsl.txt xcholsl.txx

xfredin.txx:   	xfredin.exe
		$(PRE)xfredin > xfredin.txx
		$(DIFF) xfredin.txt xfredin.txx

xmnbrak.txx:   	xmnbrak.exe
		$(PRE)xmnbrak > xmnbrak.txx
		$(DIFF) xmnbrak.txt xmnbrak.txx

xrlft3.txx:    	xrlft3.exe
		$(PRE)xrlft3 > xrlft3.txx
		$(DIFF) xrlft3.txt xrlft3.txx

xvander.txx:   	xvander.exe
		$(PRE)xvander > xvander.txx
		$(DIFF) xvander.txt xvander.txx

xchsone.txx:   	xchsone.exe
		$(PRE)xchsone > xchsone.txx
		$(DIFF) xchsone.txt xchsone.txx

xfrenel.txx:   	xfrenel.exe
		$(PRE)xfrenel > xfrenel.txx
		$(DIFF) xfrenel.txt xfrenel.txx

xmnewt.txx:    	xmnewt.exe
		$(PRE)xmnewt > xmnewt.txx
		$(DIFF) xmnewt.txt xmnewt.txx

xrofunc.txx:   	xrofunc.exe
		$(PRE)xrofunc > xrofunc.txx
		$(DIFF) xrofunc.txt xrofunc.txx

xvegas.txx:    	xvegas.exe
		$(PRE)xvegas > xvegas.txx
		$(DIFF) xvegas.txt xvegas.txx

xchstwo.txx:   	xchstwo.exe
		$(PRE)xchstwo > xchstwo.txx
		$(DIFF) xchstwo.txt xchstwo.txx

xfrprmn.txx:   	xfrprmn.exe
		$(PRE)xfrprmn > xfrprmn.txx
		$(DIFF) xfrprmn.txt xfrprmn.txx

xmoment.txx:   	xmoment.exe
		$(PRE)xmoment > xmoment.txx
		$(DIFF) xmoment.txt xmoment.txx

xrtbis.txx:    	xrtbis.exe
		$(PRE)xrtbis > xrtbis.txx
		$(DIFF) xrtbis.txt xrtbis.txx

xvoltra.txx:   	xvoltra.exe
		$(PRE)xvoltra > xvoltra.txx
		$(DIFF) xvoltra.txt xvoltra.txx

xcisi.txx:     	xcisi.exe
		$(PRE)xcisi > xcisi.txx
		$(DIFF) xcisi.txt xcisi.txx

xftest.txx:    	xftest.exe
		$(PRE)xftest > xftest.txx
		$(DIFF) xftest.txt xftest.txx

xmppi.txx:     	xmppi.exe
		$(PRE)xmppi > xmppi.txx
		$(DIFF) xmppi.txt xmppi.txx

xrtflsp.txx:   	xrtflsp.exe
		$(PRE)xrtflsp > xrtflsp.txx
		$(DIFF) xrtflsp.txt xrtflsp.txx

xwt1.txx:      	xwt1.exe
		$(PRE)xwt1 > xwt1.txx
		$(DIFF) xwt1.txt xwt1.txx

xcntab1.txx:   	xcntab1.exe
		$(PRE)xcntab1 > xcntab1.txx
		$(DIFF) xcntab1.txt xcntab1.txx

xgamdev.txx:   	xgamdev.exe
		$(PRE)xgamdev > xgamdev.txx
		$(DIFF) xgamdev.txt xgamdev.txx

xmprove.txx:   	xmprove.exe
		$(PRE)xmprove > xmprove.txx
		$(DIFF) xmprove.txt xmprove.txx

xrtnewt.txx:   	xrtnewt.exe
		$(PRE)xrtnewt > xrtnewt.txx
		$(DIFF) xrtnewt.txt xrtnewt.txx

xwtn.txx:      	xwtn.exe
		$(PRE)xwtn > xwtn.txx
		$(DIFF) xwtn.txt xwtn.txx

xcntab2.txx:   	xcntab2.exe
		$(PRE)xcntab2 > xcntab2.txx
		$(DIFF) xcntab2.txt xcntab2.txx

xgammln.txx:   	xgammln.exe
		$(PRE)xgammln > xgammln.txx
		$(DIFF) xgammln.txt xgammln.txx

xmrqcof.txx:   	xmrqcof.exe
		$(PRE)xmrqcof > xmrqcof.txx
		$(DIFF) xmrqcof.txt xmrqcof.txx

xrtsafe.txx:   	xrtsafe.exe
		$(PRE)xrtsafe > xrtsafe.txx
		$(DIFF) xrtsafe.txt xrtsafe.txx

xzbrac.txx:    	xzbrac.exe
		$(PRE)xzbrac > xzbrac.txx
		$(DIFF) xzbrac.txt xzbrac.txx

xconvlv.txx:   	xconvlv.exe
		$(PRE)xconvlv > xconvlv.txx
		$(DIFF) xconvlv.txt xconvlv.txx

xgammp.txx:    	xgammp.exe
		$(PRE)xgammp > xgammp.txx
		$(DIFF) xgammp.txt xgammp.txx

xmrqmin.txx:   	xmrqmin.exe
		$(PRE)xmrqmin > xmrqmin.txx
		$(DIFF) xmrqmin.txt xmrqmin.txx

xrtsec.txx:    	xrtsec.exe
		$(PRE)xrtsec > xrtsec.txx
		$(DIFF) xrtsec.txt xrtsec.txx

xzbrak.txx:    	xzbrak.exe
		$(PRE)xzbrak > xzbrak.txx
		$(DIFF) xzbrak.txt xzbrak.txx

xcorrel.txx:   	xcorrel.exe
		$(PRE)xcorrel > xcorrel.txx
		$(DIFF) xcorrel.txt xcorrel.txx

xgammq.txx:    	xgammq.exe
		$(PRE)xgammq > xgammq.txx
		$(DIFF) xgammq.txt xgammq.txx

xnewt.txx:     	xnewt.exe
		$(PRE)xnewt > xnewt.txx
		$(DIFF) xnewt.txt xnewt.txx

xrzextr.txx:   	xrzextr.exe
		$(PRE)xrzextr > xrzextr.txx
		$(DIFF) xrzextr.txt xrzextr.txx

xzbrent.txx:   	xzbrent.exe
		$(PRE)xzbrent > xzbrent.txx
		$(DIFF) xzbrent.txt xzbrent.txx

xcosft1.txx:   	xcosft1.exe
		$(PRE)xcosft1 > xcosft1.txx
		$(DIFF) xcosft1.txt xcosft1.txx

xgasdev.txx:   	xgasdev.exe
		$(PRE)xgasdev > xgasdev.txx
		$(DIFF) xgasdev.txt xgasdev.txx

xodeint.txx:   	xodeint.exe
		$(PRE)xodeint > xodeint.txx
		$(DIFF) xodeint.txt xodeint.txx

xsavgol.txx:   	xsavgol.exe
		$(PRE)xsavgol > xsavgol.txx
		$(DIFF) xsavgol.txt xsavgol.txx

xzrhqr.txx:    	xzrhqr.exe
		$(PRE)xzrhqr > xzrhqr.txx
		$(DIFF) xzrhqr.txt xzrhqr.txx

xcosft2.txx:   	xcosft2.exe
		$(PRE)xcosft2 > xcosft2.txx
		$(DIFF) xcosft2.txt xcosft2.txx

xgaucof.txx:   	xgaucof.exe
		$(PRE)xgaucof > xgaucof.txx
		$(DIFF) xgaucof.txt xgaucof.txx

xorthog.txx:   	xorthog.exe
		$(PRE)xorthog > xorthog.txx
		$(DIFF) xorthog.txt xorthog.txx

xscrsho.txx:   	xscrsho.exe
		$(PRE)xscrsho > xscrsho.txx
		$(DIFF) xscrsho.txt xscrsho.txx

xzriddr.txx:   	xzriddr.exe
		$(PRE)xzriddr > xzriddr.txx
		$(DIFF) xzriddr.txt xzriddr.txx

xcovsrt.txx:   	xcovsrt.exe
		$(PRE)xcovsrt > xcovsrt.txx
		$(DIFF) xcovsrt.txt xcovsrt.txx

xgauher.txx:   	xgauher.exe
		$(PRE)xgauher > xgauher.txx
		$(DIFF) xgauher.txt xgauher.txx

xpade.txx:     	xpade.exe
		$(PRE)xpade > xpade.txx
		$(DIFF) xpade.txt xpade.txx

xselect.txx:   	xselect.exe
		$(PRE)xselect > xselect.txx
		$(DIFF) xselect.txt xselect.txx

xzroots.txx:   	xzroots.exe
		$(PRE)xzroots > xzroots.txx
		$(DIFF) xzroots.txt xzroots.txx

xcrank.txx:    	xcrank.exe
		$(PRE)xcrank > xcrank.txx
		$(DIFF) xcrank.txt xcrank.txx

xgaujac.txx:   	xgaujac.exe
		$(PRE)xgaujac > xgaujac.txx
		$(DIFF) xgaujac.txt xgaujac.txx

xpccheb.txx:   	xpccheb.exe
		$(PRE)xpccheb > xpccheb.txx
		$(DIFF) xpccheb.txt xpccheb.txx

xselip.txx:    	xselip.exe
		$(PRE)xselip > xselip.txx
		$(DIFF) xselip.txt xselip.txx

