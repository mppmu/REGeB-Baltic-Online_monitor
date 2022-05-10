// :paste seg-bege-daq.scala

import scala.async.Async.{async, await}
import scala.concurrent.{Future, Promise}, scala.concurrent.duration._
import akka.actor.{Cancellable}
import daqcore.actors._, daqcore.io._, daqcore.devices._, daqcore.util._, daqcore.data._, daqcore.defaults._
import daqcore.util.fileops._

import SIS3316.dataTypes._

def exit() = { daqcoreSystem.shutdown(); daqcoreSystem.awaitTermination(); sys.exit(0); }

object logger extends Logging
logger.info("Ready")


val adc = SIS3316("vme-sis3316://minidex-fadc", "adc")

println(s"ADC identity: ${adc.identity.get}, serial number ${adc.serNo.get}")
println(s"ADC temperature: ${adc.internalTemperature.get} °C")


def configureADC_allch(): Unit = {
  adc.trigger_extern_enabled_set(allChannels --> true)
  adc.trigger_intern_enabled_set(allChannels --> false)
  adc.event_format_set(allChannels --> EventFormat())
  adc.bank_fill_threshold_stop_set(allChannels --> false)
  adc.getSync().get
}


val GCDXCore = 1
val REGeBCore = 2
val HPGeAll = Ch(GCDXCore, REGeBCore)

def configureADC_hpge(): Unit = {
  adc.trigger_intern_enabled_set(HPGeAll --> true)
  adc.input_invert_set(GCDXCore --> true)
  adc.input_invert_set(REGeBCore --> true)

  adc.trigger_gate_window_length_set(HPGeAll --> 5000)

  //adc.trigger_threshold_set(HPGeAll --> 20 * 255)
  adc.trigger_threshold_set(REGeBCore --> 20 * 255)
  adc.trigger_threshold_set(GCDXCore  --> 10 * 255)
  adc.trigger_cfd_set(HPGeAll --> CfdCtrl.CDF50Percent)
  adc.trigger_peakTime_set(HPGeAll --> 255)
  adc.trigger_gapTime_set(HPGeAll --> 250)

//  adc.energy_peakTime_set(HPGeAll --> 1500)
//  adc.energy_gapTime_set(HPGeAll --> 550)

//  adc.energy_peakTime_set(HPGeAll --> 2000)
//  adc.energy_gapTime_set(HPGeAll --> 250)

  adc.energy_peakTime_set(HPGeAll --> 1000)
  adc.energy_gapTime_set(HPGeAll --> 200)

  adc.energy_tau_table_set(HPGeAll --> 0)
  adc.energy_tau_factor_set(HPGeAll --> 34)

  adc.event_format_set(HPGeAll -->
    EventFormat(
      save_maw_values = None,
      save_energy = true,
      save_ft_maw = true,
      save_acc_78 = false,
      save_ph_acc16 = true,
      nSamples = 0,
      nMAWValues = 0
    )
  )

  adc.nsamples_pretrig_set(HPGeAll --> 3000)
  adc.nmaw_pretrig_set(HPGeAll --> 700)

  adc.bank_fill_threshold_stop_set(HPGeAll --> false)

  adc.getSync().get
  val rawEventDataSize = adc.event_format_get(HPGeAll).get vMap {_.rawEventDataSize}
  adc.bank_fill_threshold_nbytes_set(rawEventDataSize vMap {4 * _})
  adc.getSync().get
}


def configureADC(): Unit = {
  configureADC_allch()
  configureADC_hpge()
}


def printStatus() {
  println(s"Buffer count: ${adc.buffer_counter_get.get}")
  println(s"ADC temperature: ${adc.internalTemperature.get} °C")
}


configureADC()

val outputBasename="rege_bkg_starting_20170920"

adc.raw_output_file_basename_set(outputBasename)
//adc.props_output_file_basename_set(outputBasename)


def stopAfter(delay: FiniteDuration): Cancellable = scheduleOnce(delay) {
  async {
    if (await(adc.capture_enabled_get)) {
      println("Stopping ADC capture")
      adc.stopCapture()
    }
  }
}

def restartAfter(delay: FiniteDuration): Cancellable = scheduleOnce(delay) {
  async {
    if (await(adc.capture_enabled_get)) {
      println("Stopping ADC capture")
      adc.stopCapture()
      Thread.sleep(10000)
      println("Restarting ADC capture")
      adc.startCapture()
      restartAfter(delay)
    }
  }
}

def start(restartInterval: FiniteDuration = 0.seconds) = {
  adc.startCapture()
  if (restartInterval > 0.seconds) restartAfter(restartInterval)
}

def stop() = adc.stopCapture()

def runFor(time: FiniteDuration) = {
  adc.startCapture()
  stopAfter(time)
}

def runFor(time: FiniteDuration, every: FiniteDuration) = schedule(0.seconds, every) {
  adc.startCapture()
  stopAfter(time)
}


/* Usage:

start(restartInterval = 2.hours)
runFor(10.minutes)
runFor(30.seconds, 2.minutes)
runFor(2.minutes, 3.minutes)
means take data for 2 minutes, wait 1 minute (in total 3 minutes) and then take data for 2 minutes
stop()
*/


/* Filter optimization:

testTriggerMAW(adc, 1, 3000, 700, -1, -1, -1)
testTriggerMAW(adc, 1, 3000, 700, 255, 250, 2 * 255)
testTriggerMAW(adc, 1, 3000, 700, 255, 250, 400 * 255)

testEnergyMAW(adc, 1, 0, 0, -1, -1)
testEnergyMAW(adc, 1, 0, 0, 0, 39)

// gnuplot> plot "sample-values.txt" w l
// gnuplot> plot "maw-values.txt" w l
*/

