CXX = icc
CXXFLAGS = -O2 -Kc++eh -Krtti

DIFF = ./sdiff
PRE = ./
MAJOR = 101
MINOR = 0

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

everything:    	nrx_ex1 xcyclic xgaulag xpcshft xshell xairy xdawson xgauleg xpearsn xsimplx xamebsa xdbrent xgaussj xperiod xsimpr xamoeba xddpoly xgcf xpiksr2 xsinft xanneal xdecchk xgolden xpiksrt xsncndn xarcode xdf1dim xgser xplgndr xsobseq xavevar xdfpmin xhpsel xpoidev xsor xbalanc xdfridr xhpsort xpolcoe xsort xbandec xdftint xhqr xpolcof xsort2 xbanmul xdlinmin xhuffman xpoldiv xsort3 xbcucof xeclass xhunt xpolin2 xspctrm xbcuint xeclazz xhypgeo xpolint xspear xbeschb xei xicrc xpowell xsphbes xbessi xeigsrt xigray xpredic xbessi0 xelle xindexx xprobks xsplie2 xbessi1 xellf xirbit1 xpsdes xsplin2 xbessik xellpi xirbit2 xpzextr xspline xbessj xelmhes xjacobi xqgaus xsplint xbessj0 xerfcc xjulday xqrdcmp xsprsax xbessj1 xerff xkendl1 xqromb xsprsin xbessjy xerffc xkendl2 xqromo xsprspm xbessk xeulsum xks2d1s xqroot xsprstm xbessk0 xevlmem xks2d2s xqrsolv xsprstp xbessk1 xexpdev xksone xqrupdt xsprstx xbessy xexpint xkstwo xqsimp xstifbs xbessy0 xf1dim xlaguer xqtrap xstiff xbessy1 xfactln xlfit xquad3d xstoerm xbeta xfactrl xlinbcg xran xsvbksb xbetai xfasper xlinmin xran4 xsvdcmp xbico xfgauss xlocate xrank xsvdfit xbnldev xfit xlubksb xratint xsvdvar xbrent xfitexy xludcmp xratlsq xtoeplz xbroydn xfixrts xmachar xrc xtptest xbsstep xfleg xmedfit xrd xtqli xcaldat xflmoon xmemcof xrealft xtrapzd xchder xfour1 xmgfas xrf xtred2 xchebev xfourfs xmglin xrj xtridag xchebft xfourn xmidpnt xrk4 xttest xchebpc xfpoly xmiser xrkdumb xtutest xchint xfred2 xmmid xrkqs xtwofft xcholsl xfredin xmnbrak xrlft3 xvander xchsone xfrenel xmnewt xrofunc xvegas xchstwo xfrprmn xmoment xrtbis xvoltra xcisi xftest xmppi xrtflsp xwt1 xcntab1 xgamdev xmprove xrtnewt xwtn xcntab2 xgammln xmrqcof xrtsafe xzbrac xconvlv xgammp xmrqmin xrtsec xzbrak xcorrel xgammq xnewt xrzextr xzbrent xcosft1 xgasdev xodeint xsavgol xzrhqr xcosft2 xgaucof xorthog xscrsho xzriddr xcovsrt xgauher xpade xselect xzroots xcrank xgaujac xpccheb xselip 

nrx_lobj = addint.o cntab2.o frprmn.o metrop.o quad3d.o snrm.o airy.o convlv.o ftest.o mgfas.o quadct.o sobseq.o amebsa.o copy.o gamdev.o mglin.o quadmx.o solvde.o amoeba.o correl.o gammln.o midexp.o quadvl.o sor.o amotry.o cosft1.o gammp.o midinf.o ran0.o sort1.o amotsa.o cosft2.o gammq.o midpnt.o ran1.o sort2.o anneal.o covsrt.o gasdev.o midsql.o ran2.o sort3.o anorm2.o crank.o gaucof.o midsqu.o ran3.o spctrm.o arcmak.o cyclic.o gauher.o miser.o ran4.o spear.o arcode.o daub4.o gaujac.o mmid.o rank.o sphbes.o arcsum.o dawson.o gaulag.o mnbrak.o ranpt.o sphfpt.o asolve.o dbrent.o gauleg.o mnewt.o ratint.o sphoot.o atimes.o ddpoly.o gaussj.o moment.o ratlsq.o splie2.o avevar.o decchk.o gcf.o mp2dfr.o ratval.o splin2.o badluk.o df1dim.o golden.o mpdiv.o rc.o spline.o balanc.o dfpmin.o gser.o mpinv.o rd.o splint.o banbks.o dfridr.o hpsel.o mpmul.o realft.o spread.o bandec.o dftcor.o hpsort.o mpops.o rebin.o sprsax.o banmul.o dftint.o hqr.o mppi.o red.o sprsin.o bcucof.o difeq.o hufapp.o mprove.o relax.o sprspm.o bcuint.o dlinmin.o hufdec.o mpsqrt.o relax2.o sprstm.o beschb.o eclass.o hufenc.o mrqcof.o resid.o sprstp.o bessi.o eclazz.o hufmak.o mrqmin.o revcst.o sprstx.o bessi0.o ei.o hunt.o newt.o reverse.o stifbs.o bessi1.o eigsrt.o hypdrv.o odeint.o rf.o stiff.o bessik.o elle.o hypgeo.o orthog.o rj.o stoerm.o bessj.o ellf.o hypser.o pade.o rk4.o svbksb.o bessj0.o ellpi.o icrc.o pccheb.o rkck.o svdcmp.o bessj1.o elmhes.o icrc1.o pcshft.o rkdumb.o svdfit.o bessjy.o erfcc.o igray.o pearsn.o rkqs.o svdvar.o bessk.o erff.o indexx.o period.o rlft3.o toeplz.o bessk0.o erffc.o interp.o piksr2.o rofunc.o tptest.o bessk1.o eulsum.o irbit1.o piksrt.o rotate.o tqli.o bessy.o evlmem.o irbit2.o pinvs.o rsolv.o trapzd.o bessy0.o expdev.o jacobi1.o plgndr.o rstrct.o tred2.o bessy1.o expint.o jacobn.o poidev.o rtbis.o tridag.o beta.o f1dim.o julday.o polcoe.o rtflsp.o trncst.o betacf.o factln.o kendl1.o polcof.o rtnewt.o trnspt.o betai.o factrl.o kendl2.o poldiv.o rtsafe.o ttest.o bico.o fasper.o kermom.o polin2.o rtsec.o tutest.o bksub.o fdjac.o ks2d1s.o polint.o rzextr.o twofft.o bnldev.o fgauss.o ks2d2s.o powell.o savgol.o vander.o brent.o ksone.o predic.o scrsho.o vegas.o broydn.o fit.o kstwo.o probks.o select.o voltra.o bsstep.o fitexy.o laguer.o psdes.o selip.o wt1.o caldat.o fixrts.o lfit.o pwt.o sfroid.o wtn.o chder.o fleg.o linbcg.o pwtset.o shell.o wwghts.o chebev.o flmoon.o linmin.o pythag.o shoot.o zbrac.o chebft.o fmin.o lnsrch.o pzextr.o shootf.o zbrak.o chebpc.o four1.o locate.o qgaus.o simp1.o zbrent.o chint.o fourew.o lop.o qrdcmp.o simp2.o zrhqr.o chixy.o fourfs.o lubksb.o qromb.o simp3.o zriddr.o choldc.o fourn.o ludcmp.o qromo.o simplx.o zroots.o cholsl.o fpoly.o machar.o qroot.o simpr.o chsone.o fred2.o matadd.o qrsolv.o sinft.o chstwo.o fredex.o matsub.o qrupdt.o slvsm2.o cisi.o fredin.o medfit.o qsimp.o slvsml.o cntab1.o frenel.o memcof.o qtrap.o sncndn.o xsphfpt.o nruti_nm.o nrtyp_nm.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o newmatex.o bandmat.o submat.o myexcept.o

libnrx.a:      	$(nrx_lobj)
		$(AR) -cr $@ $(nrx_lobj)
		ranlib $@

nrx_ex1_obj = nrx_ex1.o cholesky.o evalue.o fft.o hholder.o jacobi.o newfft.o sort.o svd.o newmatrm.o

nrx_ex1:       	$(nrx_ex1_obj) libnrx.a
		$(CXX) -o $@ $(nrx_ex1_obj) -L. -lnrx -lm

xcyclic_obj = xcyclic.o

xcyclic:       	$(xcyclic_obj) libnrx.a
		$(CXX) -o $@ $(xcyclic_obj) -L. -lnrx -lm

xgaulag_obj = xgaulag.o

xgaulag:       	$(xgaulag_obj) libnrx.a
		$(CXX) -o $@ $(xgaulag_obj) -L. -lnrx -lm

xpcshft_obj = xpcshft.o

xpcshft:       	$(xpcshft_obj) libnrx.a
		$(CXX) -o $@ $(xpcshft_obj) -L. -lnrx -lm

xshell_obj = xshell.o

xshell:        	$(xshell_obj) libnrx.a
		$(CXX) -o $@ $(xshell_obj) -L. -lnrx -lm

xairy_obj = xairy.o

xairy:         	$(xairy_obj) libnrx.a
		$(CXX) -o $@ $(xairy_obj) -L. -lnrx -lm

xdawson_obj = xdawson.o

xdawson:       	$(xdawson_obj) libnrx.a
		$(CXX) -o $@ $(xdawson_obj) -L. -lnrx -lm

xgauleg_obj = xgauleg.o

xgauleg:       	$(xgauleg_obj) libnrx.a
		$(CXX) -o $@ $(xgauleg_obj) -L. -lnrx -lm

xpearsn_obj = xpearsn.o

xpearsn:       	$(xpearsn_obj) libnrx.a
		$(CXX) -o $@ $(xpearsn_obj) -L. -lnrx -lm

xsimplx_obj = xsimplx.o

xsimplx:       	$(xsimplx_obj) libnrx.a
		$(CXX) -o $@ $(xsimplx_obj) -L. -lnrx -lm

xamebsa_obj = xamebsa.o

xamebsa:       	$(xamebsa_obj) libnrx.a
		$(CXX) -o $@ $(xamebsa_obj) -L. -lnrx -lm

xdbrent_obj = xdbrent.o

xdbrent:       	$(xdbrent_obj) libnrx.a
		$(CXX) -o $@ $(xdbrent_obj) -L. -lnrx -lm

xgaussj_obj = xgaussj.o

xgaussj:       	$(xgaussj_obj) libnrx.a
		$(CXX) -o $@ $(xgaussj_obj) -L. -lnrx -lm

xperiod_obj = xperiod.o

xperiod:       	$(xperiod_obj) libnrx.a
		$(CXX) -o $@ $(xperiod_obj) -L. -lnrx -lm

xsimpr_obj = xsimpr.o

xsimpr:        	$(xsimpr_obj) libnrx.a
		$(CXX) -o $@ $(xsimpr_obj) -L. -lnrx -lm

xamoeba_obj = xamoeba.o

xamoeba:       	$(xamoeba_obj) libnrx.a
		$(CXX) -o $@ $(xamoeba_obj) -L. -lnrx -lm

xddpoly_obj = xddpoly.o

xddpoly:       	$(xddpoly_obj) libnrx.a
		$(CXX) -o $@ $(xddpoly_obj) -L. -lnrx -lm

xgcf_obj = xgcf.o

xgcf:          	$(xgcf_obj) libnrx.a
		$(CXX) -o $@ $(xgcf_obj) -L. -lnrx -lm

xpiksr2_obj = xpiksr2.o

xpiksr2:       	$(xpiksr2_obj) libnrx.a
		$(CXX) -o $@ $(xpiksr2_obj) -L. -lnrx -lm

xsinft_obj = xsinft.o

xsinft:        	$(xsinft_obj) libnrx.a
		$(CXX) -o $@ $(xsinft_obj) -L. -lnrx -lm

xanneal_obj = xanneal.o

xanneal:       	$(xanneal_obj) libnrx.a
		$(CXX) -o $@ $(xanneal_obj) -L. -lnrx -lm

xdecchk_obj = xdecchk.o

xdecchk:       	$(xdecchk_obj) libnrx.a
		$(CXX) -o $@ $(xdecchk_obj) -L. -lnrx -lm

xgolden_obj = xgolden.o

xgolden:       	$(xgolden_obj) libnrx.a
		$(CXX) -o $@ $(xgolden_obj) -L. -lnrx -lm

xpiksrt_obj = xpiksrt.o

xpiksrt:       	$(xpiksrt_obj) libnrx.a
		$(CXX) -o $@ $(xpiksrt_obj) -L. -lnrx -lm

xsncndn_obj = xsncndn.o

xsncndn:       	$(xsncndn_obj) libnrx.a
		$(CXX) -o $@ $(xsncndn_obj) -L. -lnrx -lm

xarcode_obj = xarcode.o

xarcode:       	$(xarcode_obj) libnrx.a
		$(CXX) -o $@ $(xarcode_obj) -L. -lnrx -lm

xdf1dim_obj = xdf1dim.o

xdf1dim:       	$(xdf1dim_obj) libnrx.a
		$(CXX) -o $@ $(xdf1dim_obj) -L. -lnrx -lm

xgser_obj = xgser.o

xgser:         	$(xgser_obj) libnrx.a
		$(CXX) -o $@ $(xgser_obj) -L. -lnrx -lm

xplgndr_obj = xplgndr.o

xplgndr:       	$(xplgndr_obj) libnrx.a
		$(CXX) -o $@ $(xplgndr_obj) -L. -lnrx -lm

xsobseq_obj = xsobseq.o

xsobseq:       	$(xsobseq_obj) libnrx.a
		$(CXX) -o $@ $(xsobseq_obj) -L. -lnrx -lm

xavevar_obj = xavevar.o

xavevar:       	$(xavevar_obj) libnrx.a
		$(CXX) -o $@ $(xavevar_obj) -L. -lnrx -lm

xdfpmin_obj = xdfpmin.o

xdfpmin:       	$(xdfpmin_obj) libnrx.a
		$(CXX) -o $@ $(xdfpmin_obj) -L. -lnrx -lm

xhpsel_obj = xhpsel.o

xhpsel:        	$(xhpsel_obj) libnrx.a
		$(CXX) -o $@ $(xhpsel_obj) -L. -lnrx -lm

xpoidev_obj = xpoidev.o

xpoidev:       	$(xpoidev_obj) libnrx.a
		$(CXX) -o $@ $(xpoidev_obj) -L. -lnrx -lm

xsor_obj = xsor.o

xsor:          	$(xsor_obj) libnrx.a
		$(CXX) -o $@ $(xsor_obj) -L. -lnrx -lm

xbalanc_obj = xbalanc.o

xbalanc:       	$(xbalanc_obj) libnrx.a
		$(CXX) -o $@ $(xbalanc_obj) -L. -lnrx -lm

xdfridr_obj = xdfridr.o

xdfridr:       	$(xdfridr_obj) libnrx.a
		$(CXX) -o $@ $(xdfridr_obj) -L. -lnrx -lm

xhpsort_obj = xhpsort.o

xhpsort:       	$(xhpsort_obj) libnrx.a
		$(CXX) -o $@ $(xhpsort_obj) -L. -lnrx -lm

xpolcoe_obj = xpolcoe.o

xpolcoe:       	$(xpolcoe_obj) libnrx.a
		$(CXX) -o $@ $(xpolcoe_obj) -L. -lnrx -lm

xsort_obj = xsort.o

xsort:         	$(xsort_obj) libnrx.a
		$(CXX) -o $@ $(xsort_obj) -L. -lnrx -lm

xbandec_obj = xbandec.o

xbandec:       	$(xbandec_obj) libnrx.a
		$(CXX) -o $@ $(xbandec_obj) -L. -lnrx -lm

xdftint_obj = xdftint.o

xdftint:       	$(xdftint_obj) libnrx.a
		$(CXX) -o $@ $(xdftint_obj) -L. -lnrx -lm

xhqr_obj = xhqr.o

xhqr:          	$(xhqr_obj) libnrx.a
		$(CXX) -o $@ $(xhqr_obj) -L. -lnrx -lm

xpolcof_obj = xpolcof.o

xpolcof:       	$(xpolcof_obj) libnrx.a
		$(CXX) -o $@ $(xpolcof_obj) -L. -lnrx -lm

xsort2_obj = xsort2.o

xsort2:        	$(xsort2_obj) libnrx.a
		$(CXX) -o $@ $(xsort2_obj) -L. -lnrx -lm

xbanmul_obj = xbanmul.o

xbanmul:       	$(xbanmul_obj) libnrx.a
		$(CXX) -o $@ $(xbanmul_obj) -L. -lnrx -lm

xdlinmin_obj = xdlinmin.o

xdlinmin:      	$(xdlinmin_obj) libnrx.a
		$(CXX) -o $@ $(xdlinmin_obj) -L. -lnrx -lm

xhuffman_obj = xhuffman.o

xhuffman:      	$(xhuffman_obj) libnrx.a
		$(CXX) -o $@ $(xhuffman_obj) -L. -lnrx -lm

xpoldiv_obj = xpoldiv.o

xpoldiv:       	$(xpoldiv_obj) libnrx.a
		$(CXX) -o $@ $(xpoldiv_obj) -L. -lnrx -lm

xsort3_obj = xsort3.o

xsort3:        	$(xsort3_obj) libnrx.a
		$(CXX) -o $@ $(xsort3_obj) -L. -lnrx -lm

xbcucof_obj = xbcucof.o

xbcucof:       	$(xbcucof_obj) libnrx.a
		$(CXX) -o $@ $(xbcucof_obj) -L. -lnrx -lm

xeclass_obj = xeclass.o

xeclass:       	$(xeclass_obj) libnrx.a
		$(CXX) -o $@ $(xeclass_obj) -L. -lnrx -lm

xhunt_obj = xhunt.o

xhunt:         	$(xhunt_obj) libnrx.a
		$(CXX) -o $@ $(xhunt_obj) -L. -lnrx -lm

xpolin2_obj = xpolin2.o

xpolin2:       	$(xpolin2_obj) libnrx.a
		$(CXX) -o $@ $(xpolin2_obj) -L. -lnrx -lm

xspctrm_obj = xspctrm.o

xspctrm:       	$(xspctrm_obj) libnrx.a
		$(CXX) -o $@ $(xspctrm_obj) -L. -lnrx -lm

xbcuint_obj = xbcuint.o

xbcuint:       	$(xbcuint_obj) libnrx.a
		$(CXX) -o $@ $(xbcuint_obj) -L. -lnrx -lm

xeclazz_obj = xeclazz.o

xeclazz:       	$(xeclazz_obj) libnrx.a
		$(CXX) -o $@ $(xeclazz_obj) -L. -lnrx -lm

xhypgeo_obj = xhypgeo.o

xhypgeo:       	$(xhypgeo_obj) libnrx.a
		$(CXX) -o $@ $(xhypgeo_obj) -L. -lnrx -lm

xpolint_obj = xpolint.o

xpolint:       	$(xpolint_obj) libnrx.a
		$(CXX) -o $@ $(xpolint_obj) -L. -lnrx -lm

xspear_obj = xspear.o

xspear:        	$(xspear_obj) libnrx.a
		$(CXX) -o $@ $(xspear_obj) -L. -lnrx -lm

xbeschb_obj = xbeschb.o

xbeschb:       	$(xbeschb_obj) libnrx.a
		$(CXX) -o $@ $(xbeschb_obj) -L. -lnrx -lm

xei_obj = xei.o

xei:           	$(xei_obj) libnrx.a
		$(CXX) -o $@ $(xei_obj) -L. -lnrx -lm

xicrc_obj = xicrc.o

xicrc:         	$(xicrc_obj) libnrx.a
		$(CXX) -o $@ $(xicrc_obj) -L. -lnrx -lm

xpowell_obj = xpowell.o

xpowell:       	$(xpowell_obj) libnrx.a
		$(CXX) -o $@ $(xpowell_obj) -L. -lnrx -lm

xsphbes_obj = xsphbes.o

xsphbes:       	$(xsphbes_obj) libnrx.a
		$(CXX) -o $@ $(xsphbes_obj) -L. -lnrx -lm

xbessi_obj = xbessi.o

xbessi:        	$(xbessi_obj) libnrx.a
		$(CXX) -o $@ $(xbessi_obj) -L. -lnrx -lm

xeigsrt_obj = xeigsrt.o

xeigsrt:       	$(xeigsrt_obj) libnrx.a
		$(CXX) -o $@ $(xeigsrt_obj) -L. -lnrx -lm

xigray_obj = xigray.o

xigray:        	$(xigray_obj) libnrx.a
		$(CXX) -o $@ $(xigray_obj) -L. -lnrx -lm

xpredic_obj = xpredic.o

xpredic:       	$(xpredic_obj) libnrx.a
		$(CXX) -o $@ $(xpredic_obj) -L. -lnrx -lm

xbessi0_obj = xbessi0.o

xbessi0:       	$(xbessi0_obj) libnrx.a
		$(CXX) -o $@ $(xbessi0_obj) -L. -lnrx -lm

xelle_obj = xelle.o

xelle:         	$(xelle_obj) libnrx.a
		$(CXX) -o $@ $(xelle_obj) -L. -lnrx -lm

xindexx_obj = xindexx.o

xindexx:       	$(xindexx_obj) libnrx.a
		$(CXX) -o $@ $(xindexx_obj) -L. -lnrx -lm

xprobks_obj = xprobks.o

xprobks:       	$(xprobks_obj) libnrx.a
		$(CXX) -o $@ $(xprobks_obj) -L. -lnrx -lm

xsplie2_obj = xsplie2.o

xsplie2:       	$(xsplie2_obj) libnrx.a
		$(CXX) -o $@ $(xsplie2_obj) -L. -lnrx -lm

xbessi1_obj = xbessi1.o

xbessi1:       	$(xbessi1_obj) libnrx.a
		$(CXX) -o $@ $(xbessi1_obj) -L. -lnrx -lm

xellf_obj = xellf.o

xellf:         	$(xellf_obj) libnrx.a
		$(CXX) -o $@ $(xellf_obj) -L. -lnrx -lm

xirbit1_obj = xirbit1.o

xirbit1:       	$(xirbit1_obj) libnrx.a
		$(CXX) -o $@ $(xirbit1_obj) -L. -lnrx -lm

xpsdes_obj = xpsdes.o

xpsdes:        	$(xpsdes_obj) libnrx.a
		$(CXX) -o $@ $(xpsdes_obj) -L. -lnrx -lm

xsplin2_obj = xsplin2.o

xsplin2:       	$(xsplin2_obj) libnrx.a
		$(CXX) -o $@ $(xsplin2_obj) -L. -lnrx -lm

xbessik_obj = xbessik.o

xbessik:       	$(xbessik_obj) libnrx.a
		$(CXX) -o $@ $(xbessik_obj) -L. -lnrx -lm

xellpi_obj = xellpi.o

xellpi:        	$(xellpi_obj) libnrx.a
		$(CXX) -o $@ $(xellpi_obj) -L. -lnrx -lm

xirbit2_obj = xirbit2.o

xirbit2:       	$(xirbit2_obj) libnrx.a
		$(CXX) -o $@ $(xirbit2_obj) -L. -lnrx -lm

xpzextr_obj = xpzextr.o

xpzextr:       	$(xpzextr_obj) libnrx.a
		$(CXX) -o $@ $(xpzextr_obj) -L. -lnrx -lm

xspline_obj = xspline.o

xspline:       	$(xspline_obj) libnrx.a
		$(CXX) -o $@ $(xspline_obj) -L. -lnrx -lm

xbessj_obj = xbessj.o

xbessj:        	$(xbessj_obj) libnrx.a
		$(CXX) -o $@ $(xbessj_obj) -L. -lnrx -lm

xelmhes_obj = xelmhes.o

xelmhes:       	$(xelmhes_obj) libnrx.a
		$(CXX) -o $@ $(xelmhes_obj) -L. -lnrx -lm

xjacobi_obj = xjacobi.o

xjacobi:       	$(xjacobi_obj) libnrx.a
		$(CXX) -o $@ $(xjacobi_obj) -L. -lnrx -lm

xqgaus_obj = xqgaus.o

xqgaus:        	$(xqgaus_obj) libnrx.a
		$(CXX) -o $@ $(xqgaus_obj) -L. -lnrx -lm

xsplint_obj = xsplint.o

xsplint:       	$(xsplint_obj) libnrx.a
		$(CXX) -o $@ $(xsplint_obj) -L. -lnrx -lm

xbessj0_obj = xbessj0.o

xbessj0:       	$(xbessj0_obj) libnrx.a
		$(CXX) -o $@ $(xbessj0_obj) -L. -lnrx -lm

xerfcc_obj = xerfcc.o

xerfcc:        	$(xerfcc_obj) libnrx.a
		$(CXX) -o $@ $(xerfcc_obj) -L. -lnrx -lm

xjulday_obj = xjulday.o

xjulday:       	$(xjulday_obj) libnrx.a
		$(CXX) -o $@ $(xjulday_obj) -L. -lnrx -lm

xqrdcmp_obj = xqrdcmp.o

xqrdcmp:       	$(xqrdcmp_obj) libnrx.a
		$(CXX) -o $@ $(xqrdcmp_obj) -L. -lnrx -lm

xsprsax_obj = xsprsax.o

xsprsax:       	$(xsprsax_obj) libnrx.a
		$(CXX) -o $@ $(xsprsax_obj) -L. -lnrx -lm

xbessj1_obj = xbessj1.o

xbessj1:       	$(xbessj1_obj) libnrx.a
		$(CXX) -o $@ $(xbessj1_obj) -L. -lnrx -lm

xerff_obj = xerff.o

xerff:         	$(xerff_obj) libnrx.a
		$(CXX) -o $@ $(xerff_obj) -L. -lnrx -lm

xkendl1_obj = xkendl1.o

xkendl1:       	$(xkendl1_obj) libnrx.a
		$(CXX) -o $@ $(xkendl1_obj) -L. -lnrx -lm

xqromb_obj = xqromb.o

xqromb:        	$(xqromb_obj) libnrx.a
		$(CXX) -o $@ $(xqromb_obj) -L. -lnrx -lm

xsprsin_obj = xsprsin.o

xsprsin:       	$(xsprsin_obj) libnrx.a
		$(CXX) -o $@ $(xsprsin_obj) -L. -lnrx -lm

xbessjy_obj = xbessjy.o

xbessjy:       	$(xbessjy_obj) libnrx.a
		$(CXX) -o $@ $(xbessjy_obj) -L. -lnrx -lm

xerffc_obj = xerffc.o

xerffc:        	$(xerffc_obj) libnrx.a
		$(CXX) -o $@ $(xerffc_obj) -L. -lnrx -lm

xkendl2_obj = xkendl2.o

xkendl2:       	$(xkendl2_obj) libnrx.a
		$(CXX) -o $@ $(xkendl2_obj) -L. -lnrx -lm

xqromo_obj = xqromo.o

xqromo:        	$(xqromo_obj) libnrx.a
		$(CXX) -o $@ $(xqromo_obj) -L. -lnrx -lm

xsprspm_obj = xsprspm.o

xsprspm:       	$(xsprspm_obj) libnrx.a
		$(CXX) -o $@ $(xsprspm_obj) -L. -lnrx -lm

xbessk_obj = xbessk.o

xbessk:        	$(xbessk_obj) libnrx.a
		$(CXX) -o $@ $(xbessk_obj) -L. -lnrx -lm

xeulsum_obj = xeulsum.o

xeulsum:       	$(xeulsum_obj) libnrx.a
		$(CXX) -o $@ $(xeulsum_obj) -L. -lnrx -lm

xks2d1s_obj = xks2d1s.o

xks2d1s:       	$(xks2d1s_obj) libnrx.a
		$(CXX) -o $@ $(xks2d1s_obj) -L. -lnrx -lm

xqroot_obj = xqroot.o

xqroot:        	$(xqroot_obj) libnrx.a
		$(CXX) -o $@ $(xqroot_obj) -L. -lnrx -lm

xsprstm_obj = xsprstm.o

xsprstm:       	$(xsprstm_obj) libnrx.a
		$(CXX) -o $@ $(xsprstm_obj) -L. -lnrx -lm

xbessk0_obj = xbessk0.o

xbessk0:       	$(xbessk0_obj) libnrx.a
		$(CXX) -o $@ $(xbessk0_obj) -L. -lnrx -lm

xevlmem_obj = xevlmem.o

xevlmem:       	$(xevlmem_obj) libnrx.a
		$(CXX) -o $@ $(xevlmem_obj) -L. -lnrx -lm

xks2d2s_obj = xks2d2s.o

xks2d2s:       	$(xks2d2s_obj) libnrx.a
		$(CXX) -o $@ $(xks2d2s_obj) -L. -lnrx -lm

xqrsolv_obj = xqrsolv.o

xqrsolv:       	$(xqrsolv_obj) libnrx.a
		$(CXX) -o $@ $(xqrsolv_obj) -L. -lnrx -lm

xsprstp_obj = xsprstp.o

xsprstp:       	$(xsprstp_obj) libnrx.a
		$(CXX) -o $@ $(xsprstp_obj) -L. -lnrx -lm

xbessk1_obj = xbessk1.o

xbessk1:       	$(xbessk1_obj) libnrx.a
		$(CXX) -o $@ $(xbessk1_obj) -L. -lnrx -lm

xexpdev_obj = xexpdev.o

xexpdev:       	$(xexpdev_obj) libnrx.a
		$(CXX) -o $@ $(xexpdev_obj) -L. -lnrx -lm

xksone_obj = xksone.o

xksone:        	$(xksone_obj) libnrx.a
		$(CXX) -o $@ $(xksone_obj) -L. -lnrx -lm

xqrupdt_obj = xqrupdt.o

xqrupdt:       	$(xqrupdt_obj) libnrx.a
		$(CXX) -o $@ $(xqrupdt_obj) -L. -lnrx -lm

xsprstx_obj = xsprstx.o

xsprstx:       	$(xsprstx_obj) libnrx.a
		$(CXX) -o $@ $(xsprstx_obj) -L. -lnrx -lm

xbessy_obj = xbessy.o

xbessy:        	$(xbessy_obj) libnrx.a
		$(CXX) -o $@ $(xbessy_obj) -L. -lnrx -lm

xexpint_obj = xexpint.o

xexpint:       	$(xexpint_obj) libnrx.a
		$(CXX) -o $@ $(xexpint_obj) -L. -lnrx -lm

xkstwo_obj = xkstwo.o

xkstwo:        	$(xkstwo_obj) libnrx.a
		$(CXX) -o $@ $(xkstwo_obj) -L. -lnrx -lm

xqsimp_obj = xqsimp.o

xqsimp:        	$(xqsimp_obj) libnrx.a
		$(CXX) -o $@ $(xqsimp_obj) -L. -lnrx -lm

xstifbs_obj = xstifbs.o

xstifbs:       	$(xstifbs_obj) libnrx.a
		$(CXX) -o $@ $(xstifbs_obj) -L. -lnrx -lm

xbessy0_obj = xbessy0.o

xbessy0:       	$(xbessy0_obj) libnrx.a
		$(CXX) -o $@ $(xbessy0_obj) -L. -lnrx -lm

xf1dim_obj = xf1dim.o

xf1dim:        	$(xf1dim_obj) libnrx.a
		$(CXX) -o $@ $(xf1dim_obj) -L. -lnrx -lm

xlaguer_obj = xlaguer.o

xlaguer:       	$(xlaguer_obj) libnrx.a
		$(CXX) -o $@ $(xlaguer_obj) -L. -lnrx -lm

xqtrap_obj = xqtrap.o

xqtrap:        	$(xqtrap_obj) libnrx.a
		$(CXX) -o $@ $(xqtrap_obj) -L. -lnrx -lm

xstiff_obj = xstiff.o

xstiff:        	$(xstiff_obj) libnrx.a
		$(CXX) -o $@ $(xstiff_obj) -L. -lnrx -lm

xbessy1_obj = xbessy1.o

xbessy1:       	$(xbessy1_obj) libnrx.a
		$(CXX) -o $@ $(xbessy1_obj) -L. -lnrx -lm

xfactln_obj = xfactln.o

xfactln:       	$(xfactln_obj) libnrx.a
		$(CXX) -o $@ $(xfactln_obj) -L. -lnrx -lm

xlfit_obj = xlfit.o

xlfit:         	$(xlfit_obj) libnrx.a
		$(CXX) -o $@ $(xlfit_obj) -L. -lnrx -lm

xquad3d_obj = xquad3d.o

xquad3d:       	$(xquad3d_obj) libnrx.a
		$(CXX) -o $@ $(xquad3d_obj) -L. -lnrx -lm

xstoerm_obj = xstoerm.o

xstoerm:       	$(xstoerm_obj) libnrx.a
		$(CXX) -o $@ $(xstoerm_obj) -L. -lnrx -lm

xbeta_obj = xbeta.o

xbeta:         	$(xbeta_obj) libnrx.a
		$(CXX) -o $@ $(xbeta_obj) -L. -lnrx -lm

xfactrl_obj = xfactrl.o

xfactrl:       	$(xfactrl_obj) libnrx.a
		$(CXX) -o $@ $(xfactrl_obj) -L. -lnrx -lm

xlinbcg_obj = xlinbcg.o

xlinbcg:       	$(xlinbcg_obj) libnrx.a
		$(CXX) -o $@ $(xlinbcg_obj) -L. -lnrx -lm

xran_obj = xran.o

xran:          	$(xran_obj) libnrx.a
		$(CXX) -o $@ $(xran_obj) -L. -lnrx -lm

xsvbksb_obj = xsvbksb.o

xsvbksb:       	$(xsvbksb_obj) libnrx.a
		$(CXX) -o $@ $(xsvbksb_obj) -L. -lnrx -lm

xbetai_obj = xbetai.o

xbetai:        	$(xbetai_obj) libnrx.a
		$(CXX) -o $@ $(xbetai_obj) -L. -lnrx -lm

xfasper_obj = xfasper.o

xfasper:       	$(xfasper_obj) libnrx.a
		$(CXX) -o $@ $(xfasper_obj) -L. -lnrx -lm

xlinmin_obj = xlinmin.o

xlinmin:       	$(xlinmin_obj) libnrx.a
		$(CXX) -o $@ $(xlinmin_obj) -L. -lnrx -lm

xran4_obj = xran4.o

xran4:         	$(xran4_obj) libnrx.a
		$(CXX) -o $@ $(xran4_obj) -L. -lnrx -lm

xsvdcmp_obj = xsvdcmp.o

xsvdcmp:       	$(xsvdcmp_obj) libnrx.a
		$(CXX) -o $@ $(xsvdcmp_obj) -L. -lnrx -lm

xbico_obj = xbico.o

xbico:         	$(xbico_obj) libnrx.a
		$(CXX) -o $@ $(xbico_obj) -L. -lnrx -lm

xfgauss_obj = xfgauss.o

xfgauss:       	$(xfgauss_obj) libnrx.a
		$(CXX) -o $@ $(xfgauss_obj) -L. -lnrx -lm

xlocate_obj = xlocate.o

xlocate:       	$(xlocate_obj) libnrx.a
		$(CXX) -o $@ $(xlocate_obj) -L. -lnrx -lm

xrank_obj = xrank.o

xrank:         	$(xrank_obj) libnrx.a
		$(CXX) -o $@ $(xrank_obj) -L. -lnrx -lm

xsvdfit_obj = xsvdfit.o

xsvdfit:       	$(xsvdfit_obj) libnrx.a
		$(CXX) -o $@ $(xsvdfit_obj) -L. -lnrx -lm

xbnldev_obj = xbnldev.o

xbnldev:       	$(xbnldev_obj) libnrx.a
		$(CXX) -o $@ $(xbnldev_obj) -L. -lnrx -lm

xfit_obj = xfit.o

xfit:          	$(xfit_obj) libnrx.a
		$(CXX) -o $@ $(xfit_obj) -L. -lnrx -lm

xlubksb_obj = xlubksb.o

xlubksb:       	$(xlubksb_obj) libnrx.a
		$(CXX) -o $@ $(xlubksb_obj) -L. -lnrx -lm

xratint_obj = xratint.o

xratint:       	$(xratint_obj) libnrx.a
		$(CXX) -o $@ $(xratint_obj) -L. -lnrx -lm

xsvdvar_obj = xsvdvar.o

xsvdvar:       	$(xsvdvar_obj) libnrx.a
		$(CXX) -o $@ $(xsvdvar_obj) -L. -lnrx -lm

xbrent_obj = xbrent.o

xbrent:        	$(xbrent_obj) libnrx.a
		$(CXX) -o $@ $(xbrent_obj) -L. -lnrx -lm

xfitexy_obj = xfitexy.o

xfitexy:       	$(xfitexy_obj) libnrx.a
		$(CXX) -o $@ $(xfitexy_obj) -L. -lnrx -lm

xludcmp_obj = xludcmp.o

xludcmp:       	$(xludcmp_obj) libnrx.a
		$(CXX) -o $@ $(xludcmp_obj) -L. -lnrx -lm

xratlsq_obj = xratlsq.o

xratlsq:       	$(xratlsq_obj) libnrx.a
		$(CXX) -o $@ $(xratlsq_obj) -L. -lnrx -lm

xtoeplz_obj = xtoeplz.o

xtoeplz:       	$(xtoeplz_obj) libnrx.a
		$(CXX) -o $@ $(xtoeplz_obj) -L. -lnrx -lm

xbroydn_obj = xbroydn.o

xbroydn:       	$(xbroydn_obj) libnrx.a
		$(CXX) -o $@ $(xbroydn_obj) -L. -lnrx -lm

xfixrts_obj = xfixrts.o

xfixrts:       	$(xfixrts_obj) libnrx.a
		$(CXX) -o $@ $(xfixrts_obj) -L. -lnrx -lm

xmachar_obj = xmachar.o

xmachar:       	$(xmachar_obj) libnrx.a
		$(CXX) -o $@ $(xmachar_obj) -L. -lnrx -lm

xrc_obj = xrc.o

xrc:           	$(xrc_obj) libnrx.a
		$(CXX) -o $@ $(xrc_obj) -L. -lnrx -lm

xtptest_obj = xtptest.o

xtptest:       	$(xtptest_obj) libnrx.a
		$(CXX) -o $@ $(xtptest_obj) -L. -lnrx -lm

xbsstep_obj = xbsstep.o

xbsstep:       	$(xbsstep_obj) libnrx.a
		$(CXX) -o $@ $(xbsstep_obj) -L. -lnrx -lm

xfleg_obj = xfleg.o

xfleg:         	$(xfleg_obj) libnrx.a
		$(CXX) -o $@ $(xfleg_obj) -L. -lnrx -lm

xmedfit_obj = xmedfit.o

xmedfit:       	$(xmedfit_obj) libnrx.a
		$(CXX) -o $@ $(xmedfit_obj) -L. -lnrx -lm

xrd_obj = xrd.o

xrd:           	$(xrd_obj) libnrx.a
		$(CXX) -o $@ $(xrd_obj) -L. -lnrx -lm

xtqli_obj = xtqli.o

xtqli:         	$(xtqli_obj) libnrx.a
		$(CXX) -o $@ $(xtqli_obj) -L. -lnrx -lm

xcaldat_obj = xcaldat.o

xcaldat:       	$(xcaldat_obj) libnrx.a
		$(CXX) -o $@ $(xcaldat_obj) -L. -lnrx -lm

xflmoon_obj = xflmoon.o

xflmoon:       	$(xflmoon_obj) libnrx.a
		$(CXX) -o $@ $(xflmoon_obj) -L. -lnrx -lm

xmemcof_obj = xmemcof.o

xmemcof:       	$(xmemcof_obj) libnrx.a
		$(CXX) -o $@ $(xmemcof_obj) -L. -lnrx -lm

xrealft_obj = xrealft.o

xrealft:       	$(xrealft_obj) libnrx.a
		$(CXX) -o $@ $(xrealft_obj) -L. -lnrx -lm

xtrapzd_obj = xtrapzd.o

xtrapzd:       	$(xtrapzd_obj) libnrx.a
		$(CXX) -o $@ $(xtrapzd_obj) -L. -lnrx -lm

xchder_obj = xchder.o

xchder:        	$(xchder_obj) libnrx.a
		$(CXX) -o $@ $(xchder_obj) -L. -lnrx -lm

xfour1_obj = xfour1.o

xfour1:        	$(xfour1_obj) libnrx.a
		$(CXX) -o $@ $(xfour1_obj) -L. -lnrx -lm

xmgfas_obj = xmgfas.o

xmgfas:        	$(xmgfas_obj) libnrx.a
		$(CXX) -o $@ $(xmgfas_obj) -L. -lnrx -lm

xrf_obj = xrf.o

xrf:           	$(xrf_obj) libnrx.a
		$(CXX) -o $@ $(xrf_obj) -L. -lnrx -lm

xtred2_obj = xtred2.o

xtred2:        	$(xtred2_obj) libnrx.a
		$(CXX) -o $@ $(xtred2_obj) -L. -lnrx -lm

xchebev_obj = xchebev.o

xchebev:       	$(xchebev_obj) libnrx.a
		$(CXX) -o $@ $(xchebev_obj) -L. -lnrx -lm

xfourfs_obj = xfourfs.o

xfourfs:       	$(xfourfs_obj) libnrx.a
		$(CXX) -o $@ $(xfourfs_obj) -L. -lnrx -lm

xmglin_obj = xmglin.o

xmglin:        	$(xmglin_obj) libnrx.a
		$(CXX) -o $@ $(xmglin_obj) -L. -lnrx -lm

xrj_obj = xrj.o

xrj:           	$(xrj_obj) libnrx.a
		$(CXX) -o $@ $(xrj_obj) -L. -lnrx -lm

xtridag_obj = xtridag.o

xtridag:       	$(xtridag_obj) libnrx.a
		$(CXX) -o $@ $(xtridag_obj) -L. -lnrx -lm

xchebft_obj = xchebft.o

xchebft:       	$(xchebft_obj) libnrx.a
		$(CXX) -o $@ $(xchebft_obj) -L. -lnrx -lm

xfourn_obj = xfourn.o

xfourn:        	$(xfourn_obj) libnrx.a
		$(CXX) -o $@ $(xfourn_obj) -L. -lnrx -lm

xmidpnt_obj = xmidpnt.o

xmidpnt:       	$(xmidpnt_obj) libnrx.a
		$(CXX) -o $@ $(xmidpnt_obj) -L. -lnrx -lm

xrk4_obj = xrk4.o

xrk4:          	$(xrk4_obj) libnrx.a
		$(CXX) -o $@ $(xrk4_obj) -L. -lnrx -lm

xttest_obj = xttest.o

xttest:        	$(xttest_obj) libnrx.a
		$(CXX) -o $@ $(xttest_obj) -L. -lnrx -lm

xchebpc_obj = xchebpc.o

xchebpc:       	$(xchebpc_obj) libnrx.a
		$(CXX) -o $@ $(xchebpc_obj) -L. -lnrx -lm

xfpoly_obj = xfpoly.o

xfpoly:        	$(xfpoly_obj) libnrx.a
		$(CXX) -o $@ $(xfpoly_obj) -L. -lnrx -lm

xmiser_obj = xmiser.o

xmiser:        	$(xmiser_obj) libnrx.a
		$(CXX) -o $@ $(xmiser_obj) -L. -lnrx -lm

xrkdumb_obj = xrkdumb.o

xrkdumb:       	$(xrkdumb_obj) libnrx.a
		$(CXX) -o $@ $(xrkdumb_obj) -L. -lnrx -lm

xtutest_obj = xtutest.o

xtutest:       	$(xtutest_obj) libnrx.a
		$(CXX) -o $@ $(xtutest_obj) -L. -lnrx -lm

xchint_obj = xchint.o

xchint:        	$(xchint_obj) libnrx.a
		$(CXX) -o $@ $(xchint_obj) -L. -lnrx -lm

xfred2_obj = xfred2.o

xfred2:        	$(xfred2_obj) libnrx.a
		$(CXX) -o $@ $(xfred2_obj) -L. -lnrx -lm

xmmid_obj = xmmid.o

xmmid:         	$(xmmid_obj) libnrx.a
		$(CXX) -o $@ $(xmmid_obj) -L. -lnrx -lm

xrkqs_obj = xrkqs.o

xrkqs:         	$(xrkqs_obj) libnrx.a
		$(CXX) -o $@ $(xrkqs_obj) -L. -lnrx -lm

xtwofft_obj = xtwofft.o

xtwofft:       	$(xtwofft_obj) libnrx.a
		$(CXX) -o $@ $(xtwofft_obj) -L. -lnrx -lm

xcholsl_obj = xcholsl.o

xcholsl:       	$(xcholsl_obj) libnrx.a
		$(CXX) -o $@ $(xcholsl_obj) -L. -lnrx -lm

xfredin_obj = xfredin.o

xfredin:       	$(xfredin_obj) libnrx.a
		$(CXX) -o $@ $(xfredin_obj) -L. -lnrx -lm

xmnbrak_obj = xmnbrak.o

xmnbrak:       	$(xmnbrak_obj) libnrx.a
		$(CXX) -o $@ $(xmnbrak_obj) -L. -lnrx -lm

xrlft3_obj = xrlft3.o

xrlft3:        	$(xrlft3_obj) libnrx.a
		$(CXX) -o $@ $(xrlft3_obj) -L. -lnrx -lm

xvander_obj = xvander.o

xvander:       	$(xvander_obj) libnrx.a
		$(CXX) -o $@ $(xvander_obj) -L. -lnrx -lm

xchsone_obj = xchsone.o

xchsone:       	$(xchsone_obj) libnrx.a
		$(CXX) -o $@ $(xchsone_obj) -L. -lnrx -lm

xfrenel_obj = xfrenel.o

xfrenel:       	$(xfrenel_obj) libnrx.a
		$(CXX) -o $@ $(xfrenel_obj) -L. -lnrx -lm

xmnewt_obj = xmnewt.o

xmnewt:        	$(xmnewt_obj) libnrx.a
		$(CXX) -o $@ $(xmnewt_obj) -L. -lnrx -lm

xrofunc_obj = xrofunc.o

xrofunc:       	$(xrofunc_obj) libnrx.a
		$(CXX) -o $@ $(xrofunc_obj) -L. -lnrx -lm

xvegas_obj = xvegas.o

xvegas:        	$(xvegas_obj) libnrx.a
		$(CXX) -o $@ $(xvegas_obj) -L. -lnrx -lm

xchstwo_obj = xchstwo.o

xchstwo:       	$(xchstwo_obj) libnrx.a
		$(CXX) -o $@ $(xchstwo_obj) -L. -lnrx -lm

xfrprmn_obj = xfrprmn.o

xfrprmn:       	$(xfrprmn_obj) libnrx.a
		$(CXX) -o $@ $(xfrprmn_obj) -L. -lnrx -lm

xmoment_obj = xmoment.o

xmoment:       	$(xmoment_obj) libnrx.a
		$(CXX) -o $@ $(xmoment_obj) -L. -lnrx -lm

xrtbis_obj = xrtbis.o

xrtbis:        	$(xrtbis_obj) libnrx.a
		$(CXX) -o $@ $(xrtbis_obj) -L. -lnrx -lm

xvoltra_obj = xvoltra.o

xvoltra:       	$(xvoltra_obj) libnrx.a
		$(CXX) -o $@ $(xvoltra_obj) -L. -lnrx -lm

xcisi_obj = xcisi.o

xcisi:         	$(xcisi_obj) libnrx.a
		$(CXX) -o $@ $(xcisi_obj) -L. -lnrx -lm

xftest_obj = xftest.o

xftest:        	$(xftest_obj) libnrx.a
		$(CXX) -o $@ $(xftest_obj) -L. -lnrx -lm

xmppi_obj = xmppi.o

xmppi:         	$(xmppi_obj) libnrx.a
		$(CXX) -o $@ $(xmppi_obj) -L. -lnrx -lm

xrtflsp_obj = xrtflsp.o

xrtflsp:       	$(xrtflsp_obj) libnrx.a
		$(CXX) -o $@ $(xrtflsp_obj) -L. -lnrx -lm

xwt1_obj = xwt1.o

xwt1:          	$(xwt1_obj) libnrx.a
		$(CXX) -o $@ $(xwt1_obj) -L. -lnrx -lm

xcntab1_obj = xcntab1.o

xcntab1:       	$(xcntab1_obj) libnrx.a
		$(CXX) -o $@ $(xcntab1_obj) -L. -lnrx -lm

xgamdev_obj = xgamdev.o

xgamdev:       	$(xgamdev_obj) libnrx.a
		$(CXX) -o $@ $(xgamdev_obj) -L. -lnrx -lm

xmprove_obj = xmprove.o

xmprove:       	$(xmprove_obj) libnrx.a
		$(CXX) -o $@ $(xmprove_obj) -L. -lnrx -lm

xrtnewt_obj = xrtnewt.o

xrtnewt:       	$(xrtnewt_obj) libnrx.a
		$(CXX) -o $@ $(xrtnewt_obj) -L. -lnrx -lm

xwtn_obj = xwtn.o

xwtn:          	$(xwtn_obj) libnrx.a
		$(CXX) -o $@ $(xwtn_obj) -L. -lnrx -lm

xcntab2_obj = xcntab2.o

xcntab2:       	$(xcntab2_obj) libnrx.a
		$(CXX) -o $@ $(xcntab2_obj) -L. -lnrx -lm

xgammln_obj = xgammln.o

xgammln:       	$(xgammln_obj) libnrx.a
		$(CXX) -o $@ $(xgammln_obj) -L. -lnrx -lm

xmrqcof_obj = xmrqcof.o

xmrqcof:       	$(xmrqcof_obj) libnrx.a
		$(CXX) -o $@ $(xmrqcof_obj) -L. -lnrx -lm

xrtsafe_obj = xrtsafe.o

xrtsafe:       	$(xrtsafe_obj) libnrx.a
		$(CXX) -o $@ $(xrtsafe_obj) -L. -lnrx -lm

xzbrac_obj = xzbrac.o

xzbrac:        	$(xzbrac_obj) libnrx.a
		$(CXX) -o $@ $(xzbrac_obj) -L. -lnrx -lm

xconvlv_obj = xconvlv.o

xconvlv:       	$(xconvlv_obj) libnrx.a
		$(CXX) -o $@ $(xconvlv_obj) -L. -lnrx -lm

xgammp_obj = xgammp.o

xgammp:        	$(xgammp_obj) libnrx.a
		$(CXX) -o $@ $(xgammp_obj) -L. -lnrx -lm

xmrqmin_obj = xmrqmin.o

xmrqmin:       	$(xmrqmin_obj) libnrx.a
		$(CXX) -o $@ $(xmrqmin_obj) -L. -lnrx -lm

xrtsec_obj = xrtsec.o

xrtsec:        	$(xrtsec_obj) libnrx.a
		$(CXX) -o $@ $(xrtsec_obj) -L. -lnrx -lm

xzbrak_obj = xzbrak.o

xzbrak:        	$(xzbrak_obj) libnrx.a
		$(CXX) -o $@ $(xzbrak_obj) -L. -lnrx -lm

xcorrel_obj = xcorrel.o

xcorrel:       	$(xcorrel_obj) libnrx.a
		$(CXX) -o $@ $(xcorrel_obj) -L. -lnrx -lm

xgammq_obj = xgammq.o

xgammq:        	$(xgammq_obj) libnrx.a
		$(CXX) -o $@ $(xgammq_obj) -L. -lnrx -lm

xnewt_obj = xnewt.o

xnewt:         	$(xnewt_obj) libnrx.a
		$(CXX) -o $@ $(xnewt_obj) -L. -lnrx -lm

xrzextr_obj = xrzextr.o

xrzextr:       	$(xrzextr_obj) libnrx.a
		$(CXX) -o $@ $(xrzextr_obj) -L. -lnrx -lm

xzbrent_obj = xzbrent.o

xzbrent:       	$(xzbrent_obj) libnrx.a
		$(CXX) -o $@ $(xzbrent_obj) -L. -lnrx -lm

xcosft1_obj = xcosft1.o

xcosft1:       	$(xcosft1_obj) libnrx.a
		$(CXX) -o $@ $(xcosft1_obj) -L. -lnrx -lm

xgasdev_obj = xgasdev.o

xgasdev:       	$(xgasdev_obj) libnrx.a
		$(CXX) -o $@ $(xgasdev_obj) -L. -lnrx -lm

xodeint_obj = xodeint.o

xodeint:       	$(xodeint_obj) libnrx.a
		$(CXX) -o $@ $(xodeint_obj) -L. -lnrx -lm

xsavgol_obj = xsavgol.o

xsavgol:       	$(xsavgol_obj) libnrx.a
		$(CXX) -o $@ $(xsavgol_obj) -L. -lnrx -lm

xzrhqr_obj = xzrhqr.o

xzrhqr:        	$(xzrhqr_obj) libnrx.a
		$(CXX) -o $@ $(xzrhqr_obj) -L. -lnrx -lm

xcosft2_obj = xcosft2.o

xcosft2:       	$(xcosft2_obj) libnrx.a
		$(CXX) -o $@ $(xcosft2_obj) -L. -lnrx -lm

xgaucof_obj = xgaucof.o

xgaucof:       	$(xgaucof_obj) libnrx.a
		$(CXX) -o $@ $(xgaucof_obj) -L. -lnrx -lm

xorthog_obj = xorthog.o

xorthog:       	$(xorthog_obj) libnrx.a
		$(CXX) -o $@ $(xorthog_obj) -L. -lnrx -lm

xscrsho_obj = xscrsho.o

xscrsho:       	$(xscrsho_obj) libnrx.a
		$(CXX) -o $@ $(xscrsho_obj) -L. -lnrx -lm

xzriddr_obj = xzriddr.o

xzriddr:       	$(xzriddr_obj) libnrx.a
		$(CXX) -o $@ $(xzriddr_obj) -L. -lnrx -lm

xcovsrt_obj = xcovsrt.o

xcovsrt:       	$(xcovsrt_obj) libnrx.a
		$(CXX) -o $@ $(xcovsrt_obj) -L. -lnrx -lm

xgauher_obj = xgauher.o

xgauher:       	$(xgauher_obj) libnrx.a
		$(CXX) -o $@ $(xgauher_obj) -L. -lnrx -lm

xpade_obj = xpade.o

xpade:         	$(xpade_obj) libnrx.a
		$(CXX) -o $@ $(xpade_obj) -L. -lnrx -lm

xselect_obj = xselect.o

xselect:       	$(xselect_obj) libnrx.a
		$(CXX) -o $@ $(xselect_obj) -L. -lnrx -lm

xzroots_obj = xzroots.o

xzroots:       	$(xzroots_obj) libnrx.a
		$(CXX) -o $@ $(xzroots_obj) -L. -lnrx -lm

xcrank_obj = xcrank.o

xcrank:        	$(xcrank_obj) libnrx.a
		$(CXX) -o $@ $(xcrank_obj) -L. -lnrx -lm

xgaujac_obj = xgaujac.o

xgaujac:       	$(xgaujac_obj) libnrx.a
		$(CXX) -o $@ $(xgaujac_obj) -L. -lnrx -lm

xpccheb_obj = xpccheb.o

xpccheb:       	$(xpccheb_obj) libnrx.a
		$(CXX) -o $@ $(xpccheb_obj) -L. -lnrx -lm

xselip_obj = xselip.o

xselip:        	$(xselip_obj) libnrx.a
		$(CXX) -o $@ $(xselip_obj) -L. -lnrx -lm

addint.o:      	addint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cntab2.o:      	cntab2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

frprmn.o:      	frprmn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

metrop.o:      	metrop.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quad3d.o:      	quad3d.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

snrm.o:        	snrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

airy.o:        	airy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

convlv.o:      	convlv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ftest.o:       	ftest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mgfas.o:       	mgfas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadct.o:      	quadct.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sobseq.o:      	sobseq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amebsa.o:      	amebsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

copy.o:        	copy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gamdev.o:      	gamdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mglin.o:       	mglin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadmx.o:      	quadmx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

solvde.o:      	solvde.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amoeba.o:      	amoeba.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

correl.o:      	correl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammln.o:      	gammln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midexp.o:      	midexp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

quadvl.o:      	quadvl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sor.o:         	sor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amotry.o:      	amotry.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cosft1.o:      	cosft1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammp.o:       	gammp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midinf.o:      	midinf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran0.o:        	ran0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort1.o:       	sort1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

amotsa.o:      	amotsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cosft2.o:      	cosft2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gammq.o:       	gammq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midpnt.o:      	midpnt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran1.o:        	ran1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort2.o:       	sort2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

anneal.o:      	anneal.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

covsrt.o:      	covsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gasdev.o:      	gasdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midsql.o:      	midsql.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran2.o:        	ran2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sort3.o:       	sort3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

anorm2.o:      	anorm2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

crank.o:       	crank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaucof.o:      	gaucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

midsqu.o:      	midsqu.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran3.o:        	ran3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spctrm.o:      	spctrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcmak.o:      	arcmak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cyclic.o:      	cyclic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gauher.o:      	gauher.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

miser.o:       	miser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ran4.o:        	ran4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spear.o:       	spear.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcode.o:      	arcode.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

daub4.o:       	daub4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaujac.o:      	gaujac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mmid.o:        	mmid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rank.o:        	rank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphbes.o:      	sphbes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

arcsum.o:      	arcsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dawson.o:      	dawson.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaulag.o:      	gaulag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mnbrak.o:      	mnbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ranpt.o:       	ranpt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphfpt.o:      	sphfpt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

asolve.o:      	asolve.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dbrent.o:      	dbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gauleg.o:      	gauleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mnewt.o:       	mnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratint.o:      	ratint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sphoot.o:      	sphoot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

atimes.o:      	atimes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ddpoly.o:      	ddpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gaussj.o:      	gaussj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

moment.o:      	moment.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratlsq.o:      	ratlsq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splie2.o:      	splie2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

avevar.o:      	avevar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

decchk.o:      	decchk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gcf.o:         	gcf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mp2dfr.o:      	mp2dfr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ratval.o:      	ratval.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splin2.o:      	splin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

badluk.o:      	badluk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

df1dim.o:      	df1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

golden.o:      	golden.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpdiv.o:       	mpdiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rc.o:          	rc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spline.o:      	spline.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

balanc.o:      	balanc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dfpmin.o:      	dfpmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

gser.o:        	gser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpinv.o:       	mpinv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rd.o:          	rd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

splint.o:      	splint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

banbks.o:      	banbks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dfridr.o:      	dfridr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hpsel.o:       	hpsel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpmul.o:       	mpmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

realft.o:      	realft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

spread.o:      	spread.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bandec.o:      	bandec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dftcor.o:      	dftcor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hpsort.o:      	hpsort.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpops.o:       	mpops.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rebin.o:       	rebin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprsax.o:      	sprsax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

banmul.o:      	banmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dftint.o:      	dftint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hqr.o:         	hqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mppi.o:        	mppi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

red.o:         	red.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprsin.o:      	sprsin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bcucof.o:      	bcucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

difeq.o:       	difeq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufapp.o:      	hufapp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mprove.o:      	mprove.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

relax.o:       	relax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprspm.o:      	sprspm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bcuint.o:      	bcuint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

dlinmin.o:     	dlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufdec.o:      	hufdec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mpsqrt.o:      	mpsqrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

relax2.o:      	relax2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstm.o:      	sprstm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

beschb.o:      	beschb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eclass.o:      	eclass.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufenc.o:      	hufenc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mrqcof.o:      	mrqcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

resid.o:       	resid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstp.o:      	sprstp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi.o:       	bessi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eclazz.o:      	eclazz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hufmak.o:      	hufmak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

mrqmin.o:      	mrqmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

revcst.o:      	revcst.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sprstx.o:      	sprstx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi0.o:      	bessi0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ei.o:          	ei.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hunt.o:        	hunt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

newt.o:        	newt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

reverse.o:     	reverse.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stifbs.o:      	stifbs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessi1.o:      	bessi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eigsrt.o:      	eigsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypdrv.o:      	hypdrv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

odeint.o:      	odeint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rf.o:          	rf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stiff.o:       	stiff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessik.o:      	bessik.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

elle.o:        	elle.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypgeo.o:      	hypgeo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

orthog.o:      	orthog.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rj.o:          	rj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

stoerm.o:      	stoerm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj.o:       	bessj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ellf.o:        	ellf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

hypser.o:      	hypser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pade.o:        	pade.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rk4.o:         	rk4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svbksb.o:      	svbksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj0.o:      	bessj0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ellpi.o:       	ellpi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

icrc.o:        	icrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pccheb.o:      	pccheb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkck.o:        	rkck.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdcmp.o:      	svdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessj1.o:      	bessj1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

elmhes.o:      	elmhes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

icrc1.o:       	icrc1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pcshft.o:      	pcshft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkdumb.o:      	rkdumb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdfit.o:      	svdfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessjy.o:      	bessjy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erfcc.o:       	erfcc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

igray.o:       	igray.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pearsn.o:      	pearsn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rkqs.o:        	rkqs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

svdvar.o:      	svdvar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk.o:       	bessk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erff.o:        	erff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

indexx.o:      	indexx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

period.o:      	period.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rlft3.o:       	rlft3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

toeplz.o:      	toeplz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk0.o:      	bessk0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

erffc.o:       	erffc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

interp.o:      	interp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

piksr2.o:      	piksr2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rofunc.o:      	rofunc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tptest.o:      	tptest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessk1.o:      	bessk1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

eulsum.o:      	eulsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

irbit1.o:      	irbit1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

piksrt.o:      	piksrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rotate.o:      	rotate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tqli.o:        	tqli.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy.o:       	bessy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

evlmem.o:      	evlmem.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

irbit2.o:      	irbit2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pinvs.o:       	pinvs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rsolv.o:       	rsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trapzd.o:      	trapzd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy0.o:      	bessy0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

expdev.o:      	expdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

jacobi1.o:     	jacobi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

plgndr.o:      	plgndr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rstrct.o:      	rstrct.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tred2.o:       	tred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bessy1.o:      	bessy1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

expint.o:      	expint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

jacobn.o:      	jacobn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

poidev.o:      	poidev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtbis.o:       	rtbis.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tridag.o:      	tridag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

beta.o:        	beta.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

f1dim.o:       	f1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

julday.o:      	julday.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polcoe.o:      	polcoe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtflsp.o:      	rtflsp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trncst.o:      	trncst.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

betacf.o:      	betacf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

factln.o:      	factln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kendl1.o:      	kendl1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polcof.o:      	polcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtnewt.o:      	rtnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

trnspt.o:      	trnspt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

betai.o:       	betai.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

factrl.o:      	factrl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kendl2.o:      	kendl2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

poldiv.o:      	poldiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtsafe.o:      	rtsafe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ttest.o:       	ttest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bico.o:        	bico.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fasper.o:      	fasper.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kermom.o:      	kermom.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polin2.o:      	polin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rtsec.o:       	rtsec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

tutest.o:      	tutest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bksub.o:       	bksub.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fdjac.o:       	fdjac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ks2d1s.o:      	ks2d1s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

polint.o:      	polint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

rzextr.o:      	rzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

twofft.o:      	twofft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bnldev.o:      	bnldev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fgauss.o:      	fgauss.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ks2d2s.o:      	ks2d2s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

powell.o:      	powell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

savgol.o:      	savgol.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

vander.o:      	vander.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

brent.o:       	brent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ksone.o:       	ksone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

predic.o:      	predic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

scrsho.o:      	scrsho.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

vegas.o:       	vegas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

broydn.o:      	broydn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fit.o:         	fit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

kstwo.o:       	kstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

probks.o:      	probks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

select.o:      	select.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

voltra.o:      	voltra.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

bsstep.o:      	bsstep.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fitexy.o:      	fitexy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

laguer.o:      	laguer.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

psdes.o:       	psdes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

selip.o:       	selip.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wt1.o:         	wt1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

caldat.o:      	caldat.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fixrts.o:      	fixrts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lfit.o:        	lfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pwt.o:         	pwt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sfroid.o:      	sfroid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wtn.o:         	wtn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chder.o:       	chder.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fleg.o:        	fleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

linbcg.o:      	linbcg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pwtset.o:      	pwtset.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shell.o:       	shell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

wwghts.o:      	wwghts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebev.o:      	chebev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

flmoon.o:      	flmoon.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

linmin.o:      	linmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pythag.o:      	pythag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shoot.o:       	shoot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrac.o:       	zbrac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebft.o:      	chebft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fmin.o:        	fmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lnsrch.o:      	lnsrch.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

pzextr.o:      	pzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

shootf.o:      	shootf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrak.o:       	zbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chebpc.o:      	chebpc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

four1.o:       	four1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

locate.o:      	locate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qgaus.o:       	qgaus.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp1.o:       	simp1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zbrent.o:      	zbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chint.o:       	chint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourew.o:      	fourew.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lop.o:         	lop.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrdcmp.o:      	qrdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp2.o:       	simp2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zrhqr.o:       	zrhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chixy.o:       	chixy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourfs.o:      	fourfs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

lubksb.o:      	lubksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qromb.o:       	qromb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simp3.o:       	simp3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zriddr.o:      	zriddr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

choldc.o:      	choldc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fourn.o:       	fourn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

ludcmp.o:      	ludcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qromo.o:       	qromo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simplx.o:      	simplx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

zroots.o:      	zroots.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cholsl.o:      	cholsl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fpoly.o:       	fpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

machar.o:      	machar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qroot.o:       	qroot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

simpr.o:       	simpr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chsone.o:      	chsone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fred2.o:       	fred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

matadd.o:      	matadd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrsolv.o:      	qrsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sinft.o:       	sinft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

chstwo.o:      	chstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fredex.o:      	fredex.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

matsub.o:      	matsub.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qrupdt.o:      	qrupdt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

slvsm2.o:      	slvsm2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cisi.o:        	cisi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

fredin.o:      	fredin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

medfit.o:      	medfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qsimp.o:       	qsimp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

slvsml.o:      	slvsml.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

cntab1.o:      	cntab1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

frenel.o:      	frenel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

memcof.o:      	memcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

qtrap.o:       	qtrap.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

sncndn.o:      	sncndn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsphfpt.o:     	xsphfpt.cpp nrtypes.h nrtyp_nm.h newmat.h nrutil.h include.h myexcept.h nruti_nm.h

nruti_nm.o:    	nruti_nm.cpp nr.h nruti_nm.h nrutil.h nrtypes.h nrtyp_nm.h newmat.h include.h myexcept.h

nrtyp_nm.o:    	nrtyp_nm.cpp nr.h nrutil_nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

newmat1.o:     	newmat1.cpp newmat.h include.h myexcept.h

newmat2.o:     	newmat2.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat3.o:     	newmat3.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat4.o:     	newmat4.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat5.o:     	newmat5.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat6.o:     	newmat6.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat7.o:     	newmat7.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

newmat8.o:     	newmat8.cpp include.h newmat.h newmatrc.h precisio.h myexcept.h controlw.h

newmatex.o:    	newmatex.cpp include.h newmat.h myexcept.h

bandmat.o:     	bandmat.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

submat.o:      	submat.cpp include.h newmat.h newmatrc.h myexcept.h controlw.h

myexcept.o:    	myexcept.cpp include.h myexcept.h

nrx_ex1.o:     	nrx_ex1.cpp newmatap.h nr.h newmat.h include.h myexcept.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h

cholesky.o:    	cholesky.cpp include.h newmat.h newmatrm.h myexcept.h

evalue.o:      	evalue.cpp include.h newmatap.h newmatrm.h precisio.h newmat.h myexcept.h

fft.o:         	fft.cpp include.h newmatap.h newmat.h myexcept.h

hholder.o:     	hholder.cpp include.h newmatap.h newmat.h myexcept.h

jacobi.o:      	jacobi.cpp include.h newmatap.h precisio.h newmatrm.h newmat.h myexcept.h

newfft.o:      	newfft.cpp newmatap.h newmat.h include.h myexcept.h

sort.o:        	sort.cpp include.h newmatap.h newmat.h myexcept.h

svd.o:         	svd.cpp include.h newmatap.h newmatrm.h precisio.h newmat.h myexcept.h

newmatrm.o:    	newmatrm.cpp newmat.h newmatrm.h include.h myexcept.h

xcyclic.o:     	xcyclic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaulag.o:     	xgaulag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpcshft.o:     	xpcshft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xshell.o:      	xshell.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xairy.o:       	xairy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdawson.o:     	xdawson.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgauleg.o:     	xgauleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpearsn.o:     	xpearsn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsimplx.o:     	xsimplx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xamebsa.o:     	xamebsa.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdbrent.o:     	xdbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaussj.o:     	xgaussj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xperiod.o:     	xperiod.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsimpr.o:      	xsimpr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xamoeba.o:     	xamoeba.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xddpoly.o:     	xddpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgcf.o:        	xgcf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpiksr2.o:     	xpiksr2.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsinft.o:      	xsinft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xanneal.o:     	xanneal.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdecchk.o:     	xdecchk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgolden.o:     	xgolden.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpiksrt.o:     	xpiksrt.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsncndn.o:     	xsncndn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xarcode.o:     	xarcode.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdf1dim.o:     	xdf1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgser.o:       	xgser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xplgndr.o:     	xplgndr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsobseq.o:     	xsobseq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xavevar.o:     	xavevar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdfpmin.o:     	xdfpmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhpsel.o:      	xhpsel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpoidev.o:     	xpoidev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsor.o:        	xsor.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbalanc.o:     	xbalanc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdfridr.o:     	xdfridr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhpsort.o:     	xhpsort.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolcoe.o:     	xpolcoe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort.o:       	xsort.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbandec.o:     	xbandec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdftint.o:     	xdftint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhqr.o:        	xhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolcof.o:     	xpolcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort2.o:      	xsort2.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbanmul.o:     	xbanmul.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xdlinmin.o:    	xdlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhuffman.o:    	xhuffman.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpoldiv.o:     	xpoldiv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsort3.o:      	xsort3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbcucof.o:     	xbcucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeclass.o:     	xeclass.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhunt.o:       	xhunt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolin2.o:     	xpolin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspctrm.o:     	xspctrm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbcuint.o:     	xbcuint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeclazz.o:     	xeclazz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xhypgeo.o:     	xhypgeo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpolint.o:     	xpolint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspear.o:      	xspear.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbeschb.o:     	xbeschb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xei.o:         	xei.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xicrc.o:       	xicrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpowell.o:     	xpowell.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsphbes.o:     	xsphbes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi.o:      	xbessi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeigsrt.o:     	xeigsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xigray.o:      	xigray.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpredic.o:     	xpredic.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi0.o:     	xbessi0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xelle.o:       	xelle.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xindexx.o:     	xindexx.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xprobks.o:     	xprobks.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplie2.o:     	xsplie2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessi1.o:     	xbessi1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xellf.o:       	xellf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xirbit1.o:     	xirbit1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpsdes.o:      	xpsdes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplin2.o:     	xsplin2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessik.o:     	xbessik.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xellpi.o:      	xellpi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xirbit2.o:     	xirbit2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpzextr.o:     	xpzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xspline.o:     	xspline.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj.o:      	xbessj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xelmhes.o:     	xelmhes.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xjacobi.o:     	xjacobi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqgaus.o:      	xqgaus.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsplint.o:     	xsplint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj0.o:     	xbessj0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerfcc.o:      	xerfcc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xjulday.o:     	xjulday.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrdcmp.o:     	xqrdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprsax.o:     	xsprsax.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessj1.o:     	xbessj1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerff.o:       	xerff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkendl1.o:     	xkendl1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqromb.o:      	xqromb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprsin.o:     	xsprsin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessjy.o:     	xbessjy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xerffc.o:      	xerffc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkendl2.o:     	xkendl2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqromo.o:      	xqromo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprspm.o:     	xsprspm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk.o:      	xbessk.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xeulsum.o:     	xeulsum.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xks2d1s.o:     	xks2d1s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqroot.o:      	xqroot.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstm.o:     	xsprstm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk0.o:     	xbessk0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xevlmem.o:     	xevlmem.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xks2d2s.o:     	xks2d2s.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrsolv.o:     	xqrsolv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstp.o:     	xsprstp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessk1.o:     	xbessk1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xexpdev.o:     	xexpdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xksone.o:      	xksone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqrupdt.o:     	xqrupdt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsprstx.o:     	xsprstx.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy.o:      	xbessy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xexpint.o:     	xexpint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xkstwo.o:      	xkstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqsimp.o:      	xqsimp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstifbs.o:     	xstifbs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy0.o:     	xbessy0.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xf1dim.o:      	xf1dim.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlaguer.o:     	xlaguer.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xqtrap.o:      	xqtrap.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstiff.o:      	xstiff.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbessy1.o:     	xbessy1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfactln.o:     	xfactln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlfit.o:       	xlfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xquad3d.o:     	xquad3d.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xstoerm.o:     	xstoerm.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbeta.o:       	xbeta.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfactrl.o:     	xfactrl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlinbcg.o:     	xlinbcg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xran.o:        	xran.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvbksb.o:     	xsvbksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbetai.o:      	xbetai.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfasper.o:     	xfasper.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlinmin.o:     	xlinmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xran4.o:       	xran4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdcmp.o:     	xsvdcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbico.o:       	xbico.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfgauss.o:     	xfgauss.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlocate.o:     	xlocate.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrank.o:       	xrank.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdfit.o:     	xsvdfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbnldev.o:     	xbnldev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfit.o:        	xfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xlubksb.o:     	xlubksb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xratint.o:     	xratint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsvdvar.o:     	xsvdvar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbrent.o:      	xbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfitexy.o:     	xfitexy.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xludcmp.o:     	xludcmp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xratlsq.o:     	xratlsq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtoeplz.o:     	xtoeplz.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbroydn.o:     	xbroydn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfixrts.o:     	xfixrts.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmachar.o:     	xmachar.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrc.o:         	xrc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtptest.o:     	xtptest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xbsstep.o:     	xbsstep.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfleg.o:       	xfleg.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmedfit.o:     	xmedfit.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrd.o:         	xrd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtqli.o:       	xtqli.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcaldat.o:     	xcaldat.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xflmoon.o:     	xflmoon.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmemcof.o:     	xmemcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrealft.o:     	xrealft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtrapzd.o:     	xtrapzd.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchder.o:      	xchder.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfour1.o:      	xfour1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmgfas.o:      	xmgfas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrf.o:         	xrf.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtred2.o:      	xtred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebev.o:     	xchebev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfourfs.o:     	xfourfs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmglin.o:      	xmglin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrj.o:         	xrj.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtridag.o:     	xtridag.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebft.o:     	xchebft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfourn.o:      	xfourn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmidpnt.o:     	xmidpnt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrk4.o:        	xrk4.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xttest.o:      	xttest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchebpc.o:     	xchebpc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfpoly.o:      	xfpoly.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmiser.o:      	xmiser.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrkdumb.o:     	xrkdumb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtutest.o:     	xtutest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchint.o:      	xchint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfred2.o:      	xfred2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmmid.o:       	xmmid.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrkqs.o:       	xrkqs.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xtwofft.o:     	xtwofft.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcholsl.o:     	xcholsl.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfredin.o:     	xfredin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmnbrak.o:     	xmnbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrlft3.o:      	xrlft3.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvander.o:     	xvander.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchsone.o:     	xchsone.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfrenel.o:     	xfrenel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmnewt.o:      	xmnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrofunc.o:     	xrofunc.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvegas.o:      	xvegas.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xchstwo.o:     	xchstwo.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xfrprmn.o:     	xfrprmn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmoment.o:     	xmoment.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtbis.o:      	xrtbis.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xvoltra.o:     	xvoltra.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcisi.o:       	xcisi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xftest.o:      	xftest.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmppi.o:       	xmppi.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtflsp.o:     	xrtflsp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xwt1.o:        	xwt1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcntab1.o:     	xcntab1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgamdev.o:     	xgamdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmprove.o:     	xmprove.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtnewt.o:     	xrtnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xwtn.o:        	xwtn.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcntab2.o:     	xcntab2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammln.o:     	xgammln.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmrqcof.o:     	xmrqcof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtsafe.o:     	xrtsafe.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrac.o:      	xzbrac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xconvlv.o:     	xconvlv.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammp.o:      	xgammp.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xmrqmin.o:     	xmrqmin.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrtsec.o:      	xrtsec.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrak.o:      	xzbrak.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcorrel.o:     	xcorrel.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgammq.o:      	xgammq.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xnewt.o:       	xnewt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xrzextr.o:     	xrzextr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzbrent.o:     	xzbrent.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcosft1.o:     	xcosft1.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgasdev.o:     	xgasdev.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xodeint.o:     	xodeint.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xsavgol.o:     	xsavgol.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzrhqr.o:      	xzrhqr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcosft2.o:     	xcosft2.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaucof.o:     	xgaucof.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xorthog.o:     	xorthog.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xscrsho.o:     	xscrsho.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzriddr.o:     	xzriddr.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcovsrt.o:     	xcovsrt.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgauher.o:     	xgauher.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpade.o:       	xpade.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xselect.o:     	xselect.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xzroots.o:     	xzroots.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xcrank.o:      	xcrank.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xgaujac.o:     	xgaujac.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xpccheb.o:     	xpccheb.cpp nr.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

xselip.o:      	xselip.cpp nr.h print_array.h nrutil.h nrtypes.h nruti_nm.h nrtyp_nm.h newmat.h include.h myexcept.h

nrx_ex1.txx:   	nrx_ex1
		$(PRE)nrx_ex1 > nrx_ex1.txx
		$(DIFF) nrx_ex1.txt nrx_ex1.txx

xcyclic.txx:   	xcyclic
		$(PRE)xcyclic > xcyclic.txx
		$(DIFF) xcyclic.txt xcyclic.txx

xgaulag.txx:   	xgaulag
		$(PRE)xgaulag > xgaulag.txx
		$(DIFF) xgaulag.txt xgaulag.txx

xpcshft.txx:   	xpcshft
		$(PRE)xpcshft > xpcshft.txx
		$(DIFF) xpcshft.txt xpcshft.txx

xshell.txx:    	xshell
		$(PRE)xshell > xshell.txx
		$(DIFF) xshell.txt xshell.txx

xairy.txx:     	xairy
		$(PRE)xairy > xairy.txx
		$(DIFF) xairy.txt xairy.txx

xdawson.txx:   	xdawson
		$(PRE)xdawson > xdawson.txx
		$(DIFF) xdawson.txt xdawson.txx

xgauleg.txx:   	xgauleg
		$(PRE)xgauleg > xgauleg.txx
		$(DIFF) xgauleg.txt xgauleg.txx

xpearsn.txx:   	xpearsn
		$(PRE)xpearsn > xpearsn.txx
		$(DIFF) xpearsn.txt xpearsn.txx

xsimplx.txx:   	xsimplx
		$(PRE)xsimplx > xsimplx.txx
		$(DIFF) xsimplx.txt xsimplx.txx

xamebsa.txx:   	xamebsa
		$(PRE)xamebsa > xamebsa.txx
		$(DIFF) xamebsa.txt xamebsa.txx

xdbrent.txx:   	xdbrent
		$(PRE)xdbrent > xdbrent.txx
		$(DIFF) xdbrent.txt xdbrent.txx

xgaussj.txx:   	xgaussj
		$(PRE)xgaussj > xgaussj.txx
		$(DIFF) xgaussj.txt xgaussj.txx

xperiod.txx:   	xperiod
		$(PRE)xperiod > xperiod.txx
		$(DIFF) xperiod.txt xperiod.txx

xsimpr.txx:    	xsimpr
		$(PRE)xsimpr > xsimpr.txx
		$(DIFF) xsimpr.txt xsimpr.txx

xamoeba.txx:   	xamoeba
		$(PRE)xamoeba > xamoeba.txx
		$(DIFF) xamoeba.txt xamoeba.txx

xddpoly.txx:   	xddpoly
		$(PRE)xddpoly > xddpoly.txx
		$(DIFF) xddpoly.txt xddpoly.txx

xgcf.txx:      	xgcf
		$(PRE)xgcf > xgcf.txx
		$(DIFF) xgcf.txt xgcf.txx

xpiksr2.txx:   	xpiksr2
		$(PRE)xpiksr2 > xpiksr2.txx
		$(DIFF) xpiksr2.txt xpiksr2.txx

xsinft.txx:    	xsinft
		$(PRE)xsinft > xsinft.txx
		$(DIFF) xsinft.txt xsinft.txx

xanneal.txx:   	xanneal
		$(PRE)xanneal > xanneal.txx
		$(DIFF) xanneal.txt xanneal.txx

xdecchk.txx:   	xdecchk
		$(PRE)xdecchk > xdecchk.txx
		$(DIFF) xdecchk.txt xdecchk.txx

xgolden.txx:   	xgolden
		$(PRE)xgolden > xgolden.txx
		$(DIFF) xgolden.txt xgolden.txx

xpiksrt.txx:   	xpiksrt
		$(PRE)xpiksrt > xpiksrt.txx
		$(DIFF) xpiksrt.txt xpiksrt.txx

xsncndn.txx:   	xsncndn
		$(PRE)xsncndn > xsncndn.txx
		$(DIFF) xsncndn.txt xsncndn.txx

xarcode.txx:   	xarcode
		$(PRE)xarcode > xarcode.txx
		$(DIFF) xarcode.txt xarcode.txx

xdf1dim.txx:   	xdf1dim
		$(PRE)xdf1dim > xdf1dim.txx
		$(DIFF) xdf1dim.txt xdf1dim.txx

xgser.txx:     	xgser
		$(PRE)xgser > xgser.txx
		$(DIFF) xgser.txt xgser.txx

xplgndr.txx:   	xplgndr
		$(PRE)xplgndr > xplgndr.txx
		$(DIFF) xplgndr.txt xplgndr.txx

xsobseq.txx:   	xsobseq
		$(PRE)xsobseq > xsobseq.txx
		$(DIFF) xsobseq.txt xsobseq.txx

xavevar.txx:   	xavevar
		$(PRE)xavevar > xavevar.txx
		$(DIFF) xavevar.txt xavevar.txx

xdfpmin.txx:   	xdfpmin
		$(PRE)xdfpmin > xdfpmin.txx
		$(DIFF) xdfpmin.txt xdfpmin.txx

xhpsel.txx:    	xhpsel
		$(PRE)xhpsel > xhpsel.txx
		$(DIFF) xhpsel.txt xhpsel.txx

xpoidev.txx:   	xpoidev
		$(PRE)xpoidev > xpoidev.txx
		$(DIFF) xpoidev.txt xpoidev.txx

xsor.txx:      	xsor
		$(PRE)xsor > xsor.txx
		$(DIFF) xsor.txt xsor.txx

xbalanc.txx:   	xbalanc
		$(PRE)xbalanc > xbalanc.txx
		$(DIFF) xbalanc.txt xbalanc.txx

xdfridr.txx:   	xdfridr
		$(PRE)xdfridr > xdfridr.txx
		$(DIFF) xdfridr.txt xdfridr.txx

xhpsort.txx:   	xhpsort
		$(PRE)xhpsort > xhpsort.txx
		$(DIFF) xhpsort.txt xhpsort.txx

xpolcoe.txx:   	xpolcoe
		$(PRE)xpolcoe > xpolcoe.txx
		$(DIFF) xpolcoe.txt xpolcoe.txx

xsort.txx:     	xsort
		$(PRE)xsort > xsort.txx
		$(DIFF) xsort.txt xsort.txx

xbandec.txx:   	xbandec
		$(PRE)xbandec > xbandec.txx
		$(DIFF) xbandec.txt xbandec.txx

xdftint.txx:   	xdftint
		$(PRE)xdftint > xdftint.txx
		$(DIFF) xdftint.txt xdftint.txx

xhqr.txx:      	xhqr
		$(PRE)xhqr > xhqr.txx
		$(DIFF) xhqr.txt xhqr.txx

xpolcof.txx:   	xpolcof
		$(PRE)xpolcof > xpolcof.txx
		$(DIFF) xpolcof.txt xpolcof.txx

xsort2.txx:    	xsort2
		$(PRE)xsort2 > xsort2.txx
		$(DIFF) xsort2.txt xsort2.txx

xbanmul.txx:   	xbanmul
		$(PRE)xbanmul > xbanmul.txx
		$(DIFF) xbanmul.txt xbanmul.txx

xdlinmin.txx:  	xdlinmin
		$(PRE)xdlinmin > xdlinmin.txx
		$(DIFF) xdlinmin.txt xdlinmin.txx

xhuffman.txx:  	xhuffman
		$(PRE)xhuffman > xhuffman.txx
		$(DIFF) xhuffman.txt xhuffman.txx

xpoldiv.txx:   	xpoldiv
		$(PRE)xpoldiv > xpoldiv.txx
		$(DIFF) xpoldiv.txt xpoldiv.txx

xsort3.txx:    	xsort3
		$(PRE)xsort3 > xsort3.txx
		$(DIFF) xsort3.txt xsort3.txx

xbcucof.txx:   	xbcucof
		$(PRE)xbcucof > xbcucof.txx
		$(DIFF) xbcucof.txt xbcucof.txx

xeclass.txx:   	xeclass
		$(PRE)xeclass > xeclass.txx
		$(DIFF) xeclass.txt xeclass.txx

xhunt.txx:     	xhunt
		$(PRE)xhunt > xhunt.txx
		$(DIFF) xhunt.txt xhunt.txx

xpolin2.txx:   	xpolin2
		$(PRE)xpolin2 > xpolin2.txx
		$(DIFF) xpolin2.txt xpolin2.txx

xspctrm.txx:   	xspctrm
		$(PRE)xspctrm > xspctrm.txx
		$(DIFF) xspctrm.txt xspctrm.txx

xbcuint.txx:   	xbcuint
		$(PRE)xbcuint > xbcuint.txx
		$(DIFF) xbcuint.txt xbcuint.txx

xeclazz.txx:   	xeclazz
		$(PRE)xeclazz > xeclazz.txx
		$(DIFF) xeclazz.txt xeclazz.txx

xhypgeo.txx:   	xhypgeo
		$(PRE)xhypgeo > xhypgeo.txx
		$(DIFF) xhypgeo.txt xhypgeo.txx

xpolint.txx:   	xpolint
		$(PRE)xpolint > xpolint.txx
		$(DIFF) xpolint.txt xpolint.txx

xspear.txx:    	xspear
		$(PRE)xspear > xspear.txx
		$(DIFF) xspear.txt xspear.txx

xbeschb.txx:   	xbeschb
		$(PRE)xbeschb > xbeschb.txx
		$(DIFF) xbeschb.txt xbeschb.txx

xei.txx:       	xei
		$(PRE)xei > xei.txx
		$(DIFF) xei.txt xei.txx

xicrc.txx:     	xicrc
		$(PRE)xicrc > xicrc.txx
		$(DIFF) xicrc.txt xicrc.txx

xpowell.txx:   	xpowell
		$(PRE)xpowell > xpowell.txx
		$(DIFF) xpowell.txt xpowell.txx

xsphbes.txx:   	xsphbes
		$(PRE)xsphbes > xsphbes.txx
		$(DIFF) xsphbes.txt xsphbes.txx

xbessi.txx:    	xbessi
		$(PRE)xbessi > xbessi.txx
		$(DIFF) xbessi.txt xbessi.txx

xeigsrt.txx:   	xeigsrt
		$(PRE)xeigsrt > xeigsrt.txx
		$(DIFF) xeigsrt.txt xeigsrt.txx

xigray.txx:    	xigray
		$(PRE)xigray > xigray.txx
		$(DIFF) xigray.txt xigray.txx

xpredic.txx:   	xpredic
		$(PRE)xpredic > xpredic.txx
		$(DIFF) xpredic.txt xpredic.txx

xbessi0.txx:   	xbessi0
		$(PRE)xbessi0 > xbessi0.txx
		$(DIFF) xbessi0.txt xbessi0.txx

xelle.txx:     	xelle
		$(PRE)xelle > xelle.txx
		$(DIFF) xelle.txt xelle.txx

xindexx.txx:   	xindexx
		$(PRE)xindexx > xindexx.txx
		$(DIFF) xindexx.txt xindexx.txx

xprobks.txx:   	xprobks
		$(PRE)xprobks > xprobks.txx
		$(DIFF) xprobks.txt xprobks.txx

xsplie2.txx:   	xsplie2
		$(PRE)xsplie2 > xsplie2.txx
		$(DIFF) xsplie2.txt xsplie2.txx

xbessi1.txx:   	xbessi1
		$(PRE)xbessi1 > xbessi1.txx
		$(DIFF) xbessi1.txt xbessi1.txx

xellf.txx:     	xellf
		$(PRE)xellf > xellf.txx
		$(DIFF) xellf.txt xellf.txx

xirbit1.txx:   	xirbit1
		$(PRE)xirbit1 > xirbit1.txx
		$(DIFF) xirbit1.txt xirbit1.txx

xpsdes.txx:    	xpsdes
		$(PRE)xpsdes > xpsdes.txx
		$(DIFF) xpsdes.txt xpsdes.txx

xsplin2.txx:   	xsplin2
		$(PRE)xsplin2 > xsplin2.txx
		$(DIFF) xsplin2.txt xsplin2.txx

xbessik.txx:   	xbessik
		$(PRE)xbessik > xbessik.txx
		$(DIFF) xbessik.txt xbessik.txx

xellpi.txx:    	xellpi
		$(PRE)xellpi > xellpi.txx
		$(DIFF) xellpi.txt xellpi.txx

xirbit2.txx:   	xirbit2
		$(PRE)xirbit2 > xirbit2.txx
		$(DIFF) xirbit2.txt xirbit2.txx

xpzextr.txx:   	xpzextr
		$(PRE)xpzextr > xpzextr.txx
		$(DIFF) xpzextr.txt xpzextr.txx

xspline.txx:   	xspline
		$(PRE)xspline > xspline.txx
		$(DIFF) xspline.txt xspline.txx

xbessj.txx:    	xbessj
		$(PRE)xbessj > xbessj.txx
		$(DIFF) xbessj.txt xbessj.txx

xelmhes.txx:   	xelmhes
		$(PRE)xelmhes > xelmhes.txx
		$(DIFF) xelmhes.txt xelmhes.txx

xjacobi.txx:   	xjacobi
		$(PRE)xjacobi > xjacobi.txx
		$(DIFF) xjacobi.txt xjacobi.txx

xqgaus.txx:    	xqgaus
		$(PRE)xqgaus > xqgaus.txx
		$(DIFF) xqgaus.txt xqgaus.txx

xsplint.txx:   	xsplint
		$(PRE)xsplint > xsplint.txx
		$(DIFF) xsplint.txt xsplint.txx

xbessj0.txx:   	xbessj0
		$(PRE)xbessj0 > xbessj0.txx
		$(DIFF) xbessj0.txt xbessj0.txx

xerfcc.txx:    	xerfcc
		$(PRE)xerfcc > xerfcc.txx
		$(DIFF) xerfcc.txt xerfcc.txx

xjulday.txx:   	xjulday
		$(PRE)xjulday > xjulday.txx
		$(DIFF) xjulday.txt xjulday.txx

xqrdcmp.txx:   	xqrdcmp
		$(PRE)xqrdcmp > xqrdcmp.txx
		$(DIFF) xqrdcmp.txt xqrdcmp.txx

xsprsax.txx:   	xsprsax
		$(PRE)xsprsax > xsprsax.txx
		$(DIFF) xsprsax.txt xsprsax.txx

xbessj1.txx:   	xbessj1
		$(PRE)xbessj1 > xbessj1.txx
		$(DIFF) xbessj1.txt xbessj1.txx

xerff.txx:     	xerff
		$(PRE)xerff > xerff.txx
		$(DIFF) xerff.txt xerff.txx

xkendl1.txx:   	xkendl1
		$(PRE)xkendl1 > xkendl1.txx
		$(DIFF) xkendl1.txt xkendl1.txx

xqromb.txx:    	xqromb
		$(PRE)xqromb > xqromb.txx
		$(DIFF) xqromb.txt xqromb.txx

xsprsin.txx:   	xsprsin
		$(PRE)xsprsin > xsprsin.txx
		$(DIFF) xsprsin.txt xsprsin.txx

xbessjy.txx:   	xbessjy
		$(PRE)xbessjy > xbessjy.txx
		$(DIFF) xbessjy.txt xbessjy.txx

xerffc.txx:    	xerffc
		$(PRE)xerffc > xerffc.txx
		$(DIFF) xerffc.txt xerffc.txx

xkendl2.txx:   	xkendl2
		$(PRE)xkendl2 > xkendl2.txx
		$(DIFF) xkendl2.txt xkendl2.txx

xqromo.txx:    	xqromo
		$(PRE)xqromo > xqromo.txx
		$(DIFF) xqromo.txt xqromo.txx

xsprspm.txx:   	xsprspm
		$(PRE)xsprspm > xsprspm.txx
		$(DIFF) xsprspm.txt xsprspm.txx

xbessk.txx:    	xbessk
		$(PRE)xbessk > xbessk.txx
		$(DIFF) xbessk.txt xbessk.txx

xeulsum.txx:   	xeulsum
		$(PRE)xeulsum > xeulsum.txx
		$(DIFF) xeulsum.txt xeulsum.txx

xks2d1s.txx:   	xks2d1s
		$(PRE)xks2d1s > xks2d1s.txx
		$(DIFF) xks2d1s.txt xks2d1s.txx

xqroot.txx:    	xqroot
		$(PRE)xqroot > xqroot.txx
		$(DIFF) xqroot.txt xqroot.txx

xsprstm.txx:   	xsprstm
		$(PRE)xsprstm > xsprstm.txx
		$(DIFF) xsprstm.txt xsprstm.txx

xbessk0.txx:   	xbessk0
		$(PRE)xbessk0 > xbessk0.txx
		$(DIFF) xbessk0.txt xbessk0.txx

xevlmem.txx:   	xevlmem
		$(PRE)xevlmem > xevlmem.txx
		$(DIFF) xevlmem.txt xevlmem.txx

xks2d2s.txx:   	xks2d2s
		$(PRE)xks2d2s > xks2d2s.txx
		$(DIFF) xks2d2s.txt xks2d2s.txx

xqrsolv.txx:   	xqrsolv
		$(PRE)xqrsolv > xqrsolv.txx
		$(DIFF) xqrsolv.txt xqrsolv.txx

xsprstp.txx:   	xsprstp
		$(PRE)xsprstp > xsprstp.txx
		$(DIFF) xsprstp.txt xsprstp.txx

xbessk1.txx:   	xbessk1
		$(PRE)xbessk1 > xbessk1.txx
		$(DIFF) xbessk1.txt xbessk1.txx

xexpdev.txx:   	xexpdev
		$(PRE)xexpdev > xexpdev.txx
		$(DIFF) xexpdev.txt xexpdev.txx

xksone.txx:    	xksone
		$(PRE)xksone > xksone.txx
		$(DIFF) xksone.txt xksone.txx

xqrupdt.txx:   	xqrupdt
		$(PRE)xqrupdt > xqrupdt.txx
		$(DIFF) xqrupdt.txt xqrupdt.txx

xsprstx.txx:   	xsprstx
		$(PRE)xsprstx > xsprstx.txx
		$(DIFF) xsprstx.txt xsprstx.txx

xbessy.txx:    	xbessy
		$(PRE)xbessy > xbessy.txx
		$(DIFF) xbessy.txt xbessy.txx

xexpint.txx:   	xexpint
		$(PRE)xexpint > xexpint.txx
		$(DIFF) xexpint.txt xexpint.txx

xkstwo.txx:    	xkstwo
		$(PRE)xkstwo > xkstwo.txx
		$(DIFF) xkstwo.txt xkstwo.txx

xqsimp.txx:    	xqsimp
		$(PRE)xqsimp > xqsimp.txx
		$(DIFF) xqsimp.txt xqsimp.txx

xstifbs.txx:   	xstifbs
		$(PRE)xstifbs > xstifbs.txx
		$(DIFF) xstifbs.txt xstifbs.txx

xbessy0.txx:   	xbessy0
		$(PRE)xbessy0 > xbessy0.txx
		$(DIFF) xbessy0.txt xbessy0.txx

xf1dim.txx:    	xf1dim
		$(PRE)xf1dim > xf1dim.txx
		$(DIFF) xf1dim.txt xf1dim.txx

xlaguer.txx:   	xlaguer
		$(PRE)xlaguer > xlaguer.txx
		$(DIFF) xlaguer.txt xlaguer.txx

xqtrap.txx:    	xqtrap
		$(PRE)xqtrap > xqtrap.txx
		$(DIFF) xqtrap.txt xqtrap.txx

xstiff.txx:    	xstiff
		$(PRE)xstiff > xstiff.txx
		$(DIFF) xstiff.txt xstiff.txx

xbessy1.txx:   	xbessy1
		$(PRE)xbessy1 > xbessy1.txx
		$(DIFF) xbessy1.txt xbessy1.txx

xfactln.txx:   	xfactln
		$(PRE)xfactln > xfactln.txx
		$(DIFF) xfactln.txt xfactln.txx

xlfit.txx:     	xlfit
		$(PRE)xlfit > xlfit.txx
		$(DIFF) xlfit.txt xlfit.txx

xquad3d.txx:   	xquad3d
		$(PRE)xquad3d > xquad3d.txx
		$(DIFF) xquad3d.txt xquad3d.txx

xstoerm.txx:   	xstoerm
		$(PRE)xstoerm > xstoerm.txx
		$(DIFF) xstoerm.txt xstoerm.txx

xbeta.txx:     	xbeta
		$(PRE)xbeta > xbeta.txx
		$(DIFF) xbeta.txt xbeta.txx

xfactrl.txx:   	xfactrl
		$(PRE)xfactrl > xfactrl.txx
		$(DIFF) xfactrl.txt xfactrl.txx

xlinbcg.txx:   	xlinbcg
		$(PRE)xlinbcg > xlinbcg.txx
		$(DIFF) xlinbcg.txt xlinbcg.txx

xran.txx:      	xran
		$(PRE)xran > xran.txx
		$(DIFF) xran.txt xran.txx

xsvbksb.txx:   	xsvbksb
		$(PRE)xsvbksb > xsvbksb.txx
		$(DIFF) xsvbksb.txt xsvbksb.txx

xbetai.txx:    	xbetai
		$(PRE)xbetai > xbetai.txx
		$(DIFF) xbetai.txt xbetai.txx

xfasper.txx:   	xfasper
		$(PRE)xfasper > xfasper.txx
		$(DIFF) xfasper.txt xfasper.txx

xlinmin.txx:   	xlinmin
		$(PRE)xlinmin > xlinmin.txx
		$(DIFF) xlinmin.txt xlinmin.txx

xran4.txx:     	xran4
		$(PRE)xran4 > xran4.txx
		$(DIFF) xran4.txt xran4.txx

xsvdcmp.txx:   	xsvdcmp
		$(PRE)xsvdcmp > xsvdcmp.txx
		$(DIFF) xsvdcmp.txt xsvdcmp.txx

xbico.txx:     	xbico
		$(PRE)xbico > xbico.txx
		$(DIFF) xbico.txt xbico.txx

xfgauss.txx:   	xfgauss
		$(PRE)xfgauss > xfgauss.txx
		$(DIFF) xfgauss.txt xfgauss.txx

xlocate.txx:   	xlocate
		$(PRE)xlocate > xlocate.txx
		$(DIFF) xlocate.txt xlocate.txx

xrank.txx:     	xrank
		$(PRE)xrank > xrank.txx
		$(DIFF) xrank.txt xrank.txx

xsvdfit.txx:   	xsvdfit
		$(PRE)xsvdfit > xsvdfit.txx
		$(DIFF) xsvdfit.txt xsvdfit.txx

xbnldev.txx:   	xbnldev
		$(PRE)xbnldev > xbnldev.txx
		$(DIFF) xbnldev.txt xbnldev.txx

xfit.txx:      	xfit
		$(PRE)xfit > xfit.txx
		$(DIFF) xfit.txt xfit.txx

xlubksb.txx:   	xlubksb
		$(PRE)xlubksb > xlubksb.txx
		$(DIFF) xlubksb.txt xlubksb.txx

xratint.txx:   	xratint
		$(PRE)xratint > xratint.txx
		$(DIFF) xratint.txt xratint.txx

xsvdvar.txx:   	xsvdvar
		$(PRE)xsvdvar > xsvdvar.txx
		$(DIFF) xsvdvar.txt xsvdvar.txx

xbrent.txx:    	xbrent
		$(PRE)xbrent > xbrent.txx
		$(DIFF) xbrent.txt xbrent.txx

xfitexy.txx:   	xfitexy
		$(PRE)xfitexy > xfitexy.txx
		$(DIFF) xfitexy.txt xfitexy.txx

xludcmp.txx:   	xludcmp
		$(PRE)xludcmp > xludcmp.txx
		$(DIFF) xludcmp.txt xludcmp.txx

xratlsq.txx:   	xratlsq
		$(PRE)xratlsq > xratlsq.txx
		$(DIFF) xratlsq.txt xratlsq.txx

xtoeplz.txx:   	xtoeplz
		$(PRE)xtoeplz > xtoeplz.txx
		$(DIFF) xtoeplz.txt xtoeplz.txx

xbroydn.txx:   	xbroydn
		$(PRE)xbroydn > xbroydn.txx
		$(DIFF) xbroydn.txt xbroydn.txx

xfixrts.txx:   	xfixrts
		$(PRE)xfixrts > xfixrts.txx
		$(DIFF) xfixrts.txt xfixrts.txx

xmachar.txx:   	xmachar
		$(PRE)xmachar > xmachar.txx
		$(DIFF) xmachar.txt xmachar.txx

xrc.txx:       	xrc
		$(PRE)xrc > xrc.txx
		$(DIFF) xrc.txt xrc.txx

xtptest.txx:   	xtptest
		$(PRE)xtptest > xtptest.txx
		$(DIFF) xtptest.txt xtptest.txx

xbsstep.txx:   	xbsstep
		$(PRE)xbsstep > xbsstep.txx
		$(DIFF) xbsstep.txt xbsstep.txx

xfleg.txx:     	xfleg
		$(PRE)xfleg > xfleg.txx
		$(DIFF) xfleg.txt xfleg.txx

xmedfit.txx:   	xmedfit
		$(PRE)xmedfit > xmedfit.txx
		$(DIFF) xmedfit.txt xmedfit.txx

xrd.txx:       	xrd
		$(PRE)xrd > xrd.txx
		$(DIFF) xrd.txt xrd.txx

xtqli.txx:     	xtqli
		$(PRE)xtqli > xtqli.txx
		$(DIFF) xtqli.txt xtqli.txx

xcaldat.txx:   	xcaldat
		$(PRE)xcaldat > xcaldat.txx
		$(DIFF) xcaldat.txt xcaldat.txx

xflmoon.txx:   	xflmoon
		$(PRE)xflmoon > xflmoon.txx
		$(DIFF) xflmoon.txt xflmoon.txx

xmemcof.txx:   	xmemcof
		$(PRE)xmemcof > xmemcof.txx
		$(DIFF) xmemcof.txt xmemcof.txx

xrealft.txx:   	xrealft
		$(PRE)xrealft > xrealft.txx
		$(DIFF) xrealft.txt xrealft.txx

xtrapzd.txx:   	xtrapzd
		$(PRE)xtrapzd > xtrapzd.txx
		$(DIFF) xtrapzd.txt xtrapzd.txx

xchder.txx:    	xchder
		$(PRE)xchder > xchder.txx
		$(DIFF) xchder.txt xchder.txx

xfour1.txx:    	xfour1
		$(PRE)xfour1 > xfour1.txx
		$(DIFF) xfour1.txt xfour1.txx

xmgfas.txx:    	xmgfas
		$(PRE)xmgfas > xmgfas.txx
		$(DIFF) xmgfas.txt xmgfas.txx

xrf.txx:       	xrf
		$(PRE)xrf > xrf.txx
		$(DIFF) xrf.txt xrf.txx

xtred2.txx:    	xtred2
		$(PRE)xtred2 > xtred2.txx
		$(DIFF) xtred2.txt xtred2.txx

xchebev.txx:   	xchebev
		$(PRE)xchebev > xchebev.txx
		$(DIFF) xchebev.txt xchebev.txx

xfourfs.txx:   	xfourfs
		$(PRE)xfourfs > xfourfs.txx
		$(DIFF) xfourfs.txt xfourfs.txx

xmglin.txx:    	xmglin
		$(PRE)xmglin > xmglin.txx
		$(DIFF) xmglin.txt xmglin.txx

xrj.txx:       	xrj
		$(PRE)xrj > xrj.txx
		$(DIFF) xrj.txt xrj.txx

xtridag.txx:   	xtridag
		$(PRE)xtridag > xtridag.txx
		$(DIFF) xtridag.txt xtridag.txx

xchebft.txx:   	xchebft
		$(PRE)xchebft > xchebft.txx
		$(DIFF) xchebft.txt xchebft.txx

xfourn.txx:    	xfourn
		$(PRE)xfourn > xfourn.txx
		$(DIFF) xfourn.txt xfourn.txx

xmidpnt.txx:   	xmidpnt
		$(PRE)xmidpnt > xmidpnt.txx
		$(DIFF) xmidpnt.txt xmidpnt.txx

xrk4.txx:      	xrk4
		$(PRE)xrk4 > xrk4.txx
		$(DIFF) xrk4.txt xrk4.txx

xttest.txx:    	xttest
		$(PRE)xttest > xttest.txx
		$(DIFF) xttest.txt xttest.txx

xchebpc.txx:   	xchebpc
		$(PRE)xchebpc > xchebpc.txx
		$(DIFF) xchebpc.txt xchebpc.txx

xfpoly.txx:    	xfpoly
		$(PRE)xfpoly > xfpoly.txx
		$(DIFF) xfpoly.txt xfpoly.txx

xmiser.txx:    	xmiser
		$(PRE)xmiser > xmiser.txx
		$(DIFF) xmiser.txt xmiser.txx

xrkdumb.txx:   	xrkdumb
		$(PRE)xrkdumb > xrkdumb.txx
		$(DIFF) xrkdumb.txt xrkdumb.txx

xtutest.txx:   	xtutest
		$(PRE)xtutest > xtutest.txx
		$(DIFF) xtutest.txt xtutest.txx

xchint.txx:    	xchint
		$(PRE)xchint > xchint.txx
		$(DIFF) xchint.txt xchint.txx

xfred2.txx:    	xfred2
		$(PRE)xfred2 > xfred2.txx
		$(DIFF) xfred2.txt xfred2.txx

xmmid.txx:     	xmmid
		$(PRE)xmmid > xmmid.txx
		$(DIFF) xmmid.txt xmmid.txx

xrkqs.txx:     	xrkqs
		$(PRE)xrkqs > xrkqs.txx
		$(DIFF) xrkqs.txt xrkqs.txx

xtwofft.txx:   	xtwofft
		$(PRE)xtwofft > xtwofft.txx
		$(DIFF) xtwofft.txt xtwofft.txx

xcholsl.txx:   	xcholsl
		$(PRE)xcholsl > xcholsl.txx
		$(DIFF) xcholsl.txt xcholsl.txx

xfredin.txx:   	xfredin
		$(PRE)xfredin > xfredin.txx
		$(DIFF) xfredin.txt xfredin.txx

xmnbrak.txx:   	xmnbrak
		$(PRE)xmnbrak > xmnbrak.txx
		$(DIFF) xmnbrak.txt xmnbrak.txx

xrlft3.txx:    	xrlft3
		$(PRE)xrlft3 > xrlft3.txx
		$(DIFF) xrlft3.txt xrlft3.txx

xvander.txx:   	xvander
		$(PRE)xvander > xvander.txx
		$(DIFF) xvander.txt xvander.txx

xchsone.txx:   	xchsone
		$(PRE)xchsone > xchsone.txx
		$(DIFF) xchsone.txt xchsone.txx

xfrenel.txx:   	xfrenel
		$(PRE)xfrenel > xfrenel.txx
		$(DIFF) xfrenel.txt xfrenel.txx

xmnewt.txx:    	xmnewt
		$(PRE)xmnewt > xmnewt.txx
		$(DIFF) xmnewt.txt xmnewt.txx

xrofunc.txx:   	xrofunc
		$(PRE)xrofunc > xrofunc.txx
		$(DIFF) xrofunc.txt xrofunc.txx

xvegas.txx:    	xvegas
		$(PRE)xvegas > xvegas.txx
		$(DIFF) xvegas.txt xvegas.txx

xchstwo.txx:   	xchstwo
		$(PRE)xchstwo > xchstwo.txx
		$(DIFF) xchstwo.txt xchstwo.txx

xfrprmn.txx:   	xfrprmn
		$(PRE)xfrprmn > xfrprmn.txx
		$(DIFF) xfrprmn.txt xfrprmn.txx

xmoment.txx:   	xmoment
		$(PRE)xmoment > xmoment.txx
		$(DIFF) xmoment.txt xmoment.txx

xrtbis.txx:    	xrtbis
		$(PRE)xrtbis > xrtbis.txx
		$(DIFF) xrtbis.txt xrtbis.txx

xvoltra.txx:   	xvoltra
		$(PRE)xvoltra > xvoltra.txx
		$(DIFF) xvoltra.txt xvoltra.txx

xcisi.txx:     	xcisi
		$(PRE)xcisi > xcisi.txx
		$(DIFF) xcisi.txt xcisi.txx

xftest.txx:    	xftest
		$(PRE)xftest > xftest.txx
		$(DIFF) xftest.txt xftest.txx

xmppi.txx:     	xmppi
		$(PRE)xmppi > xmppi.txx
		$(DIFF) xmppi.txt xmppi.txx

xrtflsp.txx:   	xrtflsp
		$(PRE)xrtflsp > xrtflsp.txx
		$(DIFF) xrtflsp.txt xrtflsp.txx

xwt1.txx:      	xwt1
		$(PRE)xwt1 > xwt1.txx
		$(DIFF) xwt1.txt xwt1.txx

xcntab1.txx:   	xcntab1
		$(PRE)xcntab1 > xcntab1.txx
		$(DIFF) xcntab1.txt xcntab1.txx

xgamdev.txx:   	xgamdev
		$(PRE)xgamdev > xgamdev.txx
		$(DIFF) xgamdev.txt xgamdev.txx

xmprove.txx:   	xmprove
		$(PRE)xmprove > xmprove.txx
		$(DIFF) xmprove.txt xmprove.txx

xrtnewt.txx:   	xrtnewt
		$(PRE)xrtnewt > xrtnewt.txx
		$(DIFF) xrtnewt.txt xrtnewt.txx

xwtn.txx:      	xwtn
		$(PRE)xwtn > xwtn.txx
		$(DIFF) xwtn.txt xwtn.txx

xcntab2.txx:   	xcntab2
		$(PRE)xcntab2 > xcntab2.txx
		$(DIFF) xcntab2.txt xcntab2.txx

xgammln.txx:   	xgammln
		$(PRE)xgammln > xgammln.txx
		$(DIFF) xgammln.txt xgammln.txx

xmrqcof.txx:   	xmrqcof
		$(PRE)xmrqcof > xmrqcof.txx
		$(DIFF) xmrqcof.txt xmrqcof.txx

xrtsafe.txx:   	xrtsafe
		$(PRE)xrtsafe > xrtsafe.txx
		$(DIFF) xrtsafe.txt xrtsafe.txx

xzbrac.txx:    	xzbrac
		$(PRE)xzbrac > xzbrac.txx
		$(DIFF) xzbrac.txt xzbrac.txx

xconvlv.txx:   	xconvlv
		$(PRE)xconvlv > xconvlv.txx
		$(DIFF) xconvlv.txt xconvlv.txx

xgammp.txx:    	xgammp
		$(PRE)xgammp > xgammp.txx
		$(DIFF) xgammp.txt xgammp.txx

xmrqmin.txx:   	xmrqmin
		$(PRE)xmrqmin > xmrqmin.txx
		$(DIFF) xmrqmin.txt xmrqmin.txx

xrtsec.txx:    	xrtsec
		$(PRE)xrtsec > xrtsec.txx
		$(DIFF) xrtsec.txt xrtsec.txx

xzbrak.txx:    	xzbrak
		$(PRE)xzbrak > xzbrak.txx
		$(DIFF) xzbrak.txt xzbrak.txx

xcorrel.txx:   	xcorrel
		$(PRE)xcorrel > xcorrel.txx
		$(DIFF) xcorrel.txt xcorrel.txx

xgammq.txx:    	xgammq
		$(PRE)xgammq > xgammq.txx
		$(DIFF) xgammq.txt xgammq.txx

xnewt.txx:     	xnewt
		$(PRE)xnewt > xnewt.txx
		$(DIFF) xnewt.txt xnewt.txx

xrzextr.txx:   	xrzextr
		$(PRE)xrzextr > xrzextr.txx
		$(DIFF) xrzextr.txt xrzextr.txx

xzbrent.txx:   	xzbrent
		$(PRE)xzbrent > xzbrent.txx
		$(DIFF) xzbrent.txt xzbrent.txx

xcosft1.txx:   	xcosft1
		$(PRE)xcosft1 > xcosft1.txx
		$(DIFF) xcosft1.txt xcosft1.txx

xgasdev.txx:   	xgasdev
		$(PRE)xgasdev > xgasdev.txx
		$(DIFF) xgasdev.txt xgasdev.txx

xodeint.txx:   	xodeint
		$(PRE)xodeint > xodeint.txx
		$(DIFF) xodeint.txt xodeint.txx

xsavgol.txx:   	xsavgol
		$(PRE)xsavgol > xsavgol.txx
		$(DIFF) xsavgol.txt xsavgol.txx

xzrhqr.txx:    	xzrhqr
		$(PRE)xzrhqr > xzrhqr.txx
		$(DIFF) xzrhqr.txt xzrhqr.txx

xcosft2.txx:   	xcosft2
		$(PRE)xcosft2 > xcosft2.txx
		$(DIFF) xcosft2.txt xcosft2.txx

xgaucof.txx:   	xgaucof
		$(PRE)xgaucof > xgaucof.txx
		$(DIFF) xgaucof.txt xgaucof.txx

xorthog.txx:   	xorthog
		$(PRE)xorthog > xorthog.txx
		$(DIFF) xorthog.txt xorthog.txx

xscrsho.txx:   	xscrsho
		$(PRE)xscrsho > xscrsho.txx
		$(DIFF) xscrsho.txt xscrsho.txx

xzriddr.txx:   	xzriddr
		$(PRE)xzriddr > xzriddr.txx
		$(DIFF) xzriddr.txt xzriddr.txx

xcovsrt.txx:   	xcovsrt
		$(PRE)xcovsrt > xcovsrt.txx
		$(DIFF) xcovsrt.txt xcovsrt.txx

xgauher.txx:   	xgauher
		$(PRE)xgauher > xgauher.txx
		$(DIFF) xgauher.txt xgauher.txx

xpade.txx:     	xpade
		$(PRE)xpade > xpade.txx
		$(DIFF) xpade.txt xpade.txx

xselect.txx:   	xselect
		$(PRE)xselect > xselect.txx
		$(DIFF) xselect.txt xselect.txx

xzroots.txx:   	xzroots
		$(PRE)xzroots > xzroots.txx
		$(DIFF) xzroots.txt xzroots.txx

xcrank.txx:    	xcrank
		$(PRE)xcrank > xcrank.txx
		$(DIFF) xcrank.txt xcrank.txx

xgaujac.txx:   	xgaujac
		$(PRE)xgaujac > xgaujac.txx
		$(DIFF) xgaujac.txt xgaujac.txx

xpccheb.txx:   	xpccheb
		$(PRE)xpccheb > xpccheb.txx
		$(DIFF) xpccheb.txt xpccheb.txx

xselip.txx:    	xselip
		$(PRE)xselip > xselip.txx
		$(DIFF) xselip.txt xselip.txx

