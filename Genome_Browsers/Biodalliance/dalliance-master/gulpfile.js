var gulp = require('gulp');
var gconcat = require('gulp-concat');
var closure = require('gulp-closure-compiler');
var browserify = require('gulp-browserify');
var rename = require('gulp-rename');
// var utils = require('gulp-util');
// var gif = require('gulp-if');

/* Want to gulp-if conditional compile everything, but currently doesn't seem to work */

gulp.task('build-worker', function() {
  gulp.src('js/fetchworker.js')
  .pipe(browserify({
    debug: true,
    nobuiltins: true
  }))
  .pipe(rename('worker-all.js'))
  .pipe(gulp.dest('build/'));
});


gulp.task('build-main', function() {
  gulp.src('js/exports.js')
  .pipe(browserify({
    debug: true,
    nobuiltins: true
  }))
  .pipe(rename('dalliance-all.js'))
  .pipe(gulp.dest('build/'));
});

gulp.task('compile-worker', function() {
  gulp.src('js/fetchworker.js')
  .pipe(browserify({
    debug: true,
    nobuiltins: true
  }))
  .pipe(rename('worker-all.js'))
  .pipe(gulp.dest('tmp/'))
  // .pipe(gif(!isDev, closure()))   // Doesn't work
  .pipe(closure())
  .pipe(gulp.dest('build/'));
});

gulp.task('compile-main', function() {
  gulp.src('js/exports.js')
  .pipe(browserify({
    debug: true,
    nobuiltins: true
  }))
  .pipe(rename('dalliance-all.js'))
  .pipe(gulp.dest('tmp/'))
  // .pipe(gif(!isDev, closure()))   // Doesn't work
  .pipe(closure())
  .pipe(gulp.dest('build/'));
});

gulp.task('watch', function() {
  gulp.watch('js/*.js', ['default']);
});

gulp.task('default', ['build-main', 'build-worker']);
gulp.task('compile', ['compile-main', 'compile-worker']);